from hashlib import new
import pandas as pd
import cobra
from cobra.io import read_sbml_model
from random import randint
import random
from corda import reaction_confidence
import re
from corda import CORDA
from cobra.flux_analysis import single_gene_deletion
import scipy.stats
import json
from cobra.flux_analysis import single_gene_deletion

from context_specific_deletions import get_all_gene_ko_reactions
from context_specific_deletions import get_gene_knockout_reactions
# State the path to the file iJO1366.xml
sbml_fname = '../models/Recon3DModel_301_no_isoforms_gene_symbols_LT15_trimed.xml'
results_path = '../results'
#Max iterations
maxit = 99 

def create_rxn_confidence(gene_confidence_dict):
    rxn_conf_dict = {}
    for r in model.reactions:
        if re.match('biomass_reaction', r.id):
            rxn_conf_dict[r.id] = 3
        else:
            if len(r.genes) > 0:
                gpr = r.gene_reaction_rule
                rxn_conf_dict[r.id] =  reaction_confidence(gpr, gene_confidence_dict)
            else:
                rxn_conf_dict[r.id] =  0
    return rxn_conf_dict


def gene_confidence(gene_dict, curr_df):
    gene_confidence_dict = {}
    confidence_zero_genes = set(model_genes) - set(gene_dict.keys())
    
    for gene in confidence_zero_genes:
        gene_confidence_dict[gene] = 0
        
    for k,v in gene_dict.items():
        gene_expr_stats_df = curr_df.loc[k]
        q_1 = gene_expr_stats_df['p1']
        q_2 = gene_expr_stats_df['p2']
        q_4 = gene_expr_stats_df['p3']
        if not v == 'nan':
            if v <= q_1:
                gene_confidence_dict[k] = -1
            elif v <= q_2:
                gene_confidence_dict[k] = 1
            elif v < q_4:
                gene_confidence_dict[k] = 2
            elif v >= q_4:
                gene_confidence_dict[k] = 3
        else:
            #If 'nan' then assign confidence == 0
            gene_confidence_dict[k] = 0
    
    return gene_confidence_dict


def curr_df(x,y,z):
    aux_df = pd.DataFrame()
    aux_df['p1'] = df_expression.quantile(x/100)
    aux_df['p2'] = df_expression.quantile(y/100)
    aux_df['p3'] = df_expression.quantile(z/100)
    aux_df = aux_df.fillna(0)
    print(aux_df.head())
    return aux_df

def knockouts(model):
    knockout = single_gene_deletion(model, gene_list=model_genes)

# This is a fix to get the gene's id as the index
    knockout['ids'] = [list(i)[0] for i in knockout.ids]
    knockout = knockout.set_index('ids')

    return knockout

def run_context_specific_ko(model, conf_genes_dict, conf_threshold=2,
                            objectives=('biomass_reaction')):
    result_dict = {}
    model_genes = {g.id for g in model.genes}
    objectives = [r for r in model.reactions if r.id in objectives]
    for r in model.reactions:
        if r.objective_coefficient == 0:
            continue
        r.objective_coefficient = 0

    for obj_rxn in objectives:
        result_dict[obj_rxn.id] = {}
        obj_rxn.objective_coefficient = 1
        solution = model.optimize()
        result_dict[obj_rxn.id]['wild_type'] = solution.objective_value
        obj_rxn.objective_coefficient = 0

    gene_ko_reactions = get_gene_knockout_reactions(model)

    all_genes_ko_reactions = get_all_gene_ko_reactions(model, conf_genes_dict,
                                                       conf_threshold=conf_threshold)

    # Set of genes identified by contextualized gpr evaluation
    cs_gene_ko = set(all_genes_ko_reactions.keys()) - set(gene_ko_reactions.keys())
    for gene, rxn_list in all_genes_ko_reactions.items():
        bounds_dict = {}
        for rxn in rxn_list:
            bounds_dict[rxn.id] = rxn.bounds
            rxn.bounds = (0, 0)

        for obj_rxn in objectives:
            obj_rxn.objective_coefficient = 1
            solution = model.optimize()
            result_dict[obj_rxn.id][gene] = solution.objective_value
            obj_rxn.objective_coefficient = 0

        for rxn in rxn_list:
            rxn.bounds = bounds_dict[rxn.id]

    result_dict['confidence'] = {}
    result_dict['in_model'] = {}
    result_dict['is_cs_ko'] = {}
    result_dict['inactivate_reactions'] = {}
    genes_inactivate_reactions = set(all_genes_ko_reactions.keys())
    for gene, conf in conf_genes_dict.items():
        result_dict['confidence'][gene] = conf
        result_dict['is_cs_ko'][gene] = False
        if gene in genes_inactivate_reactions:
            result_dict['inactivate_reactions'][gene] = len(all_genes_ko_reactions[gene])
            if gene in cs_gene_ko:
                result_dict['is_cs_ko'][gene] = True
        else:
            result_dict['inactivate_reactions'][gene] = 0
            result_dict['biomass_reaction'][gene] = result_dict['biomass_reaction']['wild_type']
            #result_dict['ATPM'][gene] = result_dict['ATPM']['wild_type']

        if gene in model_genes:
            result_dict['in_model'][gene] = True
        else:
            result_dict['in_model'][gene] = False

    df_results = pd.DataFrame(result_dict)
    df_results.index.name = 'gene_id'

    return df_results

def pval(x,y,z):
    aux_dict = {}
    print('Creating current df:', end=' ')
    aux_df = curr_df(x,y,z)
    print('OK')
    print('Creating gene confidence dict:', end=' ')
    gene_confidence_dict = gene_confidence(gene_dict, aux_df)
    with open(f'../models/ga_models_2/gene_conf_{x}_{y}_{z}.json', 'w') as fp:
        json.dump(gene_confidence_dict, fp)
    print('OK')
    print('Creating rxn confidence dict:', end=' ')
    rxn_conf_dict = create_rxn_confidence(gene_confidence_dict)
    with open(f'../models/ga_models_2/rxn_conf_{x}_{y}_{z}.json', 'w') as fp:
        json.dump(rxn_conf_dict, fp)
    print('OK')
    print('Building corda model:', end=' ')
    opt = CORDA(model, rxn_conf_dict)
    opt.build()
    print(opt)
    print('OK')
    model_cs = opt.cobra_model("model_build")
    model_cs.objective = "biomass_reaction"
    cobra.io.write_sbml_model(model_cs, f'../models/ga_models_2/Recon3D_LT15_trimed_corda_model_{x}_{y}_{z}.xml')
    #knockout = knockouts(model_cs)
    ko_df = run_context_specific_ko(model, gene_confidence_dict)
    knockout = ko_df['biomass_reaction']
    print(knockout.head())
    solution_cs = model_cs.optimize()
    objective_value_cs = solution_cs.objective_value
    ko_df.to_parquet(f'../models/ga_models_2/kos_{x}_{y}_{z}.parquet')
    # We define a threshold to define whether the reduction on the biomass flux is considered lethal.
    threshold = objective_value_cs*0.01

    # Use this threshold to find which set of genes' knock out reduce the predicted growth below the threshold.
    insilico_lethals_cs = set(knockout.index[knockout< threshold])
    # The set of non-essential genes are the genes with a growth value above the threshold.
    insilico_non_lethals_cs = set(knockout.index[knockout> threshold])
    lethal_ceres_cs_df = ceres_stats_df.loc[ceres_stats_df.index.intersection(insilico_lethals_cs)]
    non_lethal_ceres_cs_df = ceres_stats_df.loc[ceres_stats_df.index.intersection(insilico_non_lethals_cs)] 


    result = scipy.stats.mannwhitneyu(lethal_ceres_cs_df['mean_ceres'],non_lethal_ceres_cs_df['mean_ceres'])
    aux_dict['pval'] = result.pvalue
    aux_dict['essentials'] = len(insilico_lethals_cs)
    aux_dict['non_essentials'] = len(insilico_non_lethals_cs)
    with open(f'../models/ga_models_2/info_{x}_{y}_{z}.json', 'w') as f:
        json.dump(aux_dict, f)

    return result.pvalue


def fitness(x,y,z):
    ans = pval(x,y,z)
    if ans == 0:
        return 99999
    else:
        return abs(1/ans)

# Read the model
print('Reading the model', end=' ')
global model
model = read_sbml_model(sbml_fname)
model.solver = 'glpk'
model.objective = "biomass_reaction"
solution = model.optimize()
objective_value = solution.objective_value
print('Done')

cell_line = "SW620_LARGE_INTESTINE"
global model_genes
model_genes = [g.id for g in model.genes]

print('Reading dataframes', end=' ')
expression_ceres_df = pd.read_parquet('../data/processed/expression_ceres_df.parquet')
expression_stats_df = pd.read_parquet('../data/processed/expression_stats_df.parquet')
global ceres_stats_df
ceres_stats_df = pd.read_parquet('../data/processed/ceres_stats_df.parquet')
global df_expression
df_expression = pd.read_parquet('../data/processed/CCLE_expression_processed.parquet')
print('OK')

global gene_dict
gene_dict = expression_ceres_df.to_dict()['expression']


#Generate solutions
print('Generating random solutions:', end=' ')
solutions = []
for s in range(99):
    solutions.append((randint(5,25),randint(26,74), randint(75,95)))
print('OK')

for i in range(maxit):
    rankedsolutions = []
    for s in solutions:
        rankedsolutions.append( (fitness(s[0], s[1], s[2]), s) )
    rankedsolutions.sort()
    rankedsolutions.reverse()

    print(f'=== Gen {i} best solutions ===')
    print(rankedsolutions[0])

    bestsolutions = rankedsolutions[:10]

    elements = []
    for e in bestsolutions:
        elements.append(e[1][0])
        elements.append(e[1][1])
        elements.append(e[1][2])

    newGen= []
    for _ in range(99):
        e1 = random.choice(elements) + randint(-1,1)
        e2 = random.choice(elements) + randint(-1,1)
        e3 = random.choice(elements) + randint(-1,1)

        newGen.append((e1,e2,e3))
    
    solutions = newGen
