import pandas as pd
import cobra
from cobra.io import read_sbml_model
import re
from cobra.flux_analysis import single_gene_deletion
import os
# State the path to the file iJO1366.xml
sbml_fname = '../models/'

# Read the models
for file in os.listdir(sbml_fname):
    if re.match('Recon3D_LT15_trimed_corda_model_',file):
        print(file)
        model = read_sbml_model(f'{sbml_fname}/{file}')
        model.objective = "biomass_reaction"
        model.solver = 'glpk'
        solution = model.optimize()
        objective_value = solution.objective_value

        # Then get the list of all the genes
        all_genes = [g.id for g in model.genes]

        # Running in silico (takes a while)
        knockout = single_gene_deletion(model, gene_list=all_genes)

        # This is a fix to get the gene's id as the index
        knockout['ids'] = [list(i)[0] for i in knockout.ids]
        knockout = knockout.set_index('ids')

        # The output of the function single_gene_deletion is a dataframe
        knockout.sort_values(by='growth')
        # We define a threshold to define whether the reduction on the biomass flux is considered lethal.
        threshold = objective_value*0.01

        # Use this threshold to find which set of genes' knock out reduce the predicted growth below the threshold.
        insilico_lethals = set(knockout.index[knockout.growth < threshold])
        # The set of non-essential genes are the genes with a growth value above the threshold.
        insilico_non_lethals = set(knockout.index[knockout.growth > threshold])

        print("in-silico lethals:", len(insilico_lethals))
        print("in-silico non lethals:", len(insilico_non_lethals))