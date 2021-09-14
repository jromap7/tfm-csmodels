import numpy as np
import pandas as pd
from cobra.io import read_sbml_model
from cobra.util import solvers
import json
from cobra.flux_analysis import single_gene_deletion

from context_specific_deletions import get_all_gene_ko_reactions
from context_specific_deletions import get_gene_knockout_reactions



######################PARAMETERS######################################
#Path of processed data
processed_path = './data/processed/experimental_cs_models'
#Name of teh Recon3D model to import
sbml_fname = './models/Recon3DModel_301_no_isoforms_gene_symbols_rpmi_trimed.xml'
#Name of the cell line we will work with
cell_line = 'SW620_LARGE_INTESTINE'

results = './results'
######################################################################


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


treatments = ['25_50_75_ccle']
for treat in treatments:
    fname = f'{processed_path}/{treat}_model.xml'
    model = read_sbml_model(fname)
    model.solver = 'glpk'

    gene_fname = f'{processed_path}/{treat}_gene_conf.json'
    with open(gene_fname) as f:
        gene_conf_dict = json.load(f)

    rxn_fname = f'{treat}_rxn_conf.json'
    with open(f'{processed_path}/{rxn_fname}') as f:
        rxn_conf_dict = json.load(f)
    print(f'Computing {treat} model...', end='')
    df = run_context_specific_ko(model, gene_conf_dict)
    #df = single_gene_deletion(model)
    print('Done.')
    print(f'Saving results from {treat} model into folder: {results}...')
    df.to_parquet(f'{results}/{treat}_single_gene_del.parquet')
    print('Done.')
