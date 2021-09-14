import pandas as pd
import re, json
import cobra
from cobra.io import read_sbml_model, write_sbml_model
from corda import CORDA
from corda import reaction_confidence

print('Reading Recon3D model rpmi trimed...', end='')
sbml_fname = '../models/Recon3DModel_301_no_isoforms_gene_symbols_rpmi_trimed.xml'

# Read the model
model = read_sbml_model(sbml_fname)
model.solver = 'glpk'
model.objective = "biomass_reaction"
model_genes = [g.id for g in model.genes]

print('Done')

print('Reading Experimental array data...', end='')

df = pd.read_parquet('../data/processed/GSE56496_geneid.parquet')
df = df.rename(columns={'GSM1362700':'control1', 'GSM1362701':'30ug_1', 'GSM1362702':'60ug_1',
'GSM1362703':'100ug_1', 'GSM1362704':'control2', 'GSM1362705':'30ug_2', 'GSM1362706':'60ug_2', 'GSM1362707':'100ug_2'})
df = df.loc[df.index.intersection(model_genes)]
df['control'] = (df['control1'] + df['control2']) / 2
df['30ug'] = (df['30ug_1'] + df['30ug_2']) / 2
df['60ug'] = (df['60ug_1'] + df['60ug_2']) / 2
df['100ug'] = (df['100ug_1'] + df['100ug_2']) / 2

print('Done')

print('Processing percentile expression tresholds...', end='')

df_expression = pd.read_parquet('../data/processed/CCLE_expression_processed.parquet')
aux_df = pd.DataFrame()
aux_df['p8'] = df_expression.quantile(8/100)
aux_df['p60'] = df_expression.quantile(60/100)
aux_df['p95'] = df_expression.quantile(95/100)
aux_df = aux_df.fillna(0)
aux_df = aux_df.loc[aux_df.index.intersection(model_genes)]
aux_df = aux_df.apply(lambda x: x*1.37 + 6.64)

ccle_df = pd.DataFrame()
ccle_df['mean'] = df_expression.mean()
ccle_df['p8'] = df_expression.quantile(8/100)
ccle_df['p60'] = df_expression.quantile(60/100)
ccle_df['p95'] = df_expression.quantile(95/100)
ccle_df = ccle_df.fillna(0)
ccle_df = ccle_df.loc[ccle_df.index.intersection(model_genes)]

common_genes = list(set(df.index) & set(aux_df.index))
aux_df = aux_df.loc[aux_df.index.intersection(common_genes)]
df = df.loc[df.index.intersection(common_genes)]
print('Done')

def gene_confidence(gene_dict, curr_df):
    gene_confidence_dict = {}
    confidence_zero_genes = set(model_genes) - set(gene_dict.keys())
    
    for gene in confidence_zero_genes:
        gene_confidence_dict[gene] = 0
        
    for k,v in gene_dict.items():
        print(curr_df.head())
        gene_expr_stats_df = curr_df.loc[k]
        q_1 = gene_expr_stats_df['p8']
        q_2 = gene_expr_stats_df['p60']
        q_4 = gene_expr_stats_df['p95']
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

#treatment = ['30ug', '60ug', '100ug']
treatment = ['ccle']
for treat in treatment:
    if treat == 'ccle':
        gene_dict = ccle_df.to_dict()['mean']
    else:
        gene_dict = df.to_dict()[treat]
    print('Computing gene confidences...', end='')
    if treat == 'ccle':
        gene_confidence_dict = gene_confidence(gene_dict, ccle_df)
    else:    
        gene_confidence_dict = gene_confidence(gene_dict, aux_df)
    print('Done')
    print(f'Saving {treat} gene confidence as {treat}_gene_conf.json:', end='')
    with open(f'../data/processed/experimental_cs_models/{treat}_gene_conf.json', 'w') as fp:
        json.dump(gene_confidence_dict, fp)
    print('Done.')
    print('Computing rxn confidences...', end='')
    rxn_conf_dict = create_rxn_confidence(gene_confidence_dict)
    print('Done')
    print(f'Saving {treat} rxn confidence as {treat}_gene_conf.json:', end='')
    with open(f'../data/processed/experimental_cs_models/{treat}_rxn_conf.json', 'w') as fp:
        json.dump(rxn_conf_dict, fp)
    print('Done.')
    print(f'Building corda for the {treat} model:', end=' ')
    opt = CORDA(model, rxn_conf_dict)
    opt.build()
    print(opt)
    print('OK')
    model_cs = opt.cobra_model("model_build")
    model_cs.objective = "biomass_reaction"
    write_sbml_model(model_cs, f'../data/processed/experimental_cs_models/{treat}_cs_model.xml')