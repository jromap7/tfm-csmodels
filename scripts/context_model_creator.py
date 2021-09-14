import pandas as pd
import cobra
from cobra.io import read_sbml_model
from corda import reaction_confidence, CORDA
from cobra.flux_analysis import single_gene_deletion
import re

# Fun to generate gene_confidence dictionary from gene_dict of model (and expression)
def gene_confidence(gene_dict, test_vec):
    gene_confidence_dict = {}
    global model_genes
    for k,v in gene_dict.items():
        gene_expr_stats_df = expression_stats_df.loc[k]
        q_1 = gene_expr_stats_df[test_vec[0]]
        q_2 = gene_expr_stats_df[test_vec[1]]
        q_3 = gene_expr_stats_df[test_vec[2]]
        q_4 = gene_expr_stats_df[test_vec[3]]
        if not v == 'nan':
            if v < q_1:
                gene_confidence_dict[k] = -1
            elif v < q_1:
                gene_confidence_dict[k] = 1
            elif v < q_2:
                gene_confidence_dict[k] = 2
            elif v >= q_3:
                gene_confidence_dict[k] = 3
        else:
            #If 'nan' then assign confidence == 0
            gene_confidence_dict[k] = 0
    
    return gene_confidence_dict

# State the path to the file iJO1366.xml
sbml_fname = '../models/Recon3DModel_301_no_isoforms_gene_symbols_LT15_trimed.xml'

# Read the model
model = read_sbml_model(sbml_fname)

model.objective = "biomass_reaction"

model.solver = 'glpk'
solution = model.optimize()
objective_value = solution.objective_value

expression_ceres_df = pd.read_parquet('../data/processed/expression_ceres_df.parquet')
expression_stats_df = pd.read_parquet('../data/processed/expression_stats_df.parquet')
ceres_stats_df = pd.read_parquet('../data/processed/ceres_stats_df.parquet')

# Then get the list of all the genes
all_genes = [g.id for g in model.genes]

# Running in silico (takes a while)
knockout = single_gene_deletion(model, gene_list=all_genes)

# This is a fix to get the gene's id as the index
knockout['ids'] = [list(i)[0] for i in knockout.ids]
knockout = knockout.set_index('ids')


# We define a threshold to define whether the reduction on the biomass flux is considered lethal.
threshold = objective_value*0.999

# Use this threshold to find which set of genes' knock out reduce the predicted growth below the threshold.
insilico_lethals = set(knockout.index[knockout.growth < threshold])
# The set of non-essential genes are the genes with a growth value above the threshold.
insilico_non_lethals = set(knockout.index[knockout.growth > threshold])

print("in-silico lethals:", len(insilico_lethals))
print("in-silico non lethals:", len(insilico_non_lethals))


lethal_ceres_df = ceres_stats_df.loc[ceres_stats_df.index.intersection(insilico_lethals)]
non_lethal_ceres_df = ceres_stats_df.loc[ceres_stats_df.index.intersection(insilico_non_lethals)]

gene_dict = expression_ceres_df.to_dict()['expression']

vect = [['q10_exp', 'q25_expr', 'q50_expr', 'q75_expr'],['q25_expr', 'q50_expr', 'q75_expr', 'q90_expr'], ['q10_expr', 'q50_expr', 'q75_expr', 'q90_expr'], ['q10_expr', 'q25_expr', 'q75_expr', 'q90_expr'], ['q10_expr', 'q25_expr', 'q50_expr', 'q90_expr']]

for test_vect in vect:
    gene_confidence_dict = gene_confidence(gene_dict, test_vect)
    rxn_conf_dict = {}
    for r in model.reactions:
        if re.match('biomass', r.id):
            rxn_conf_dict[r.id] = 3
        else:
            if len(r.genes) > 0:
                gpr = r.gene_reaction_rule
                rxn_conf_dict[r.id] =  reaction_confidence(gpr, gene_confidence_dict)
            else:
                rxn_conf_dict[r.id] =  0

    opt = CORDA(model, rxn_conf_dict)
    opt.build()
    print(opt)
    model_cs = opt.cobra_model("model_build")
    model_cs.objective = "biomass_reaction"
    cobra.io.write_sbml_model(model_cs, f'../models/Recon3D_LT15_trimed_corda_model_{test_vect[0]}_{test_vect[1]}_{test_vect[2]}_{test_vect[3]}.xml')