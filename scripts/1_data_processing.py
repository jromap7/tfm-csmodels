import cobra
from cobra.io import read_sbml_model
import pandas as pd
import re

######################PARAMETERS######################################
#Path of processed data
processed_path = './data/processed'
#Raw data path
raw_path = './data/raw'
#Name of teh Recon3D model to import
sbml_fname = './models/Recon3DModel_301_no_isoforms_gene_symbols_rpmi_trimed.xml'
#Name of the cell line we will work with
cell_line = 'SW620_LARGE_INTESTINE'
######################################################################

# Read the model
print('Reading the Recon3D Model...')
model = read_sbml_model(sbml_fname)
model.solver = 'glpk'
model_genes = set([gene.id for gene in model.genes])
#print(model_genes)
print('Done.')
# Read the sample info file and load it as pandas DF
print('Reading sample_info.csv...')
df_sample_info = pd.read_csv(f'{raw_path}/sample_info.csv', index_col=0)

print(df_sample_info.head())
ccle_names = df_sample_info.to_dict()['CCLE_Name']
df_sample_info.to_parquet(f'{processed_path}/sample_info.parquet')
print('Done.')
#print(ccle_names)

# Read the CCLE expression file and load it as pandas DF
print('Reading CCLE_expression.csv...')
df_expression = pd.read_csv(f'{raw_path}/CCLE_expression.csv', index_col=0)


#Rename index with appropriate names instead of IDs
df_expression = df_expression.rename(ccle_names, axis=0)
#Rename column index without the ID numbers; only leaving gene IDs
column_rename = {elem: re.sub(r' \([0-9]+\)$',r'', elem) for elem in df_expression.columns}
df_expression = df_expression.rename(column_rename, axis=1)
df_expression = df_expression.fillna(0)
#df_expression = df_expression.loc[df_expression.index.intersection(model_genes)]
print('Done.')

print(f'Saving modified CCLE_expression_processed into {processed_path}')
df_expression.to_parquet(f'{processed_path}/CCLE_expression_processed.parquet')
#print(df_expression.head())
print('Done.')

print('Comparing gene IDs of the dataset wih gene IDs present in the Recon3D model...')
exp_genes = set(df_expression.columns.values)
#print(exp_genes)
print('Finding gene IDs present in both dataset and model...')
intersect_genes = list(model_genes & exp_genes)
print('Done.')
#print(intersect_genes)

#Processing Achilles
achilles_df = pd.read_csv('./data/raw/Achilles_gene_effect.csv', index_col=0)
#Rename the DepMap_ID
achilles_df = achilles_df.rename(ccle_names, axis=0)
#Rename the genes (substract the brackets)
achilles_df = achilles_df.rename({elem:elem.split(' ')[0] for elem in achilles_df.columns}, axis=1)
achilles_df.to_parquet(f'{processed_path}/achilles_df.parquet')

#Creating dataframe of the cell line with gene expression and ceres
pd1 = pd.DataFrame(df_expression.loc[cell_line])
pd1 = pd1.rename({'SW620_LARGE_INTESTINE':'expression'}, axis=1)
pd2 = pd.DataFrame(achilles_df.loc[cell_line])
pd2 = pd2.rename({'SW620_LARGE_INTESTINE':'ceres'}, axis=1)
expression_ceres_df = pd.merge(pd1, pd2, left_index=True, right_index=True, how='outer')
expression_ceres_df = expression_ceres_df.loc[expression_ceres_df.index.intersection(intersect_genes)]
expression_ceres_df = expression_ceres_df.fillna(0)


expression_ceres_df.to_parquet(f'{processed_path}/expression_ceres_df.parquet')

expression_stats_df = pd.DataFrame()

expression_stats_df['mean_expression'] = df_expression.mean(axis=0)
expression_stats_df['median_expression'] = df_expression.median(axis=0)
expression_stats_df['q10_expr'] = df_expression.quantile(0.1)
expression_stats_df['q25_expr'] = df_expression.quantile(0.25)
expression_stats_df['q50_expr'] = df_expression.quantile(0.5)
expression_stats_df['q75_expr'] = df_expression.quantile(0.75)
expression_stats_df['q90_expr'] = df_expression.quantile(0.9)
expression_stats_df['min'] = df_expression.min()
expression_stats_df['max'] = df_expression.max()
expression_stats_df = expression_stats_df.loc[expression_stats_df.index.intersection(intersect_genes)]
expression_stats_df.to_parquet(f'{processed_path}/expression_stats_df.parquet')

ceres_stats_df = pd.DataFrame()

ceres_stats_df['mean_ceres'] = achilles_df.mean(axis=0)
ceres_stats_df['median_ceres'] = achilles_df.median(axis=0)
ceres_stats_df['q10_ceres'] = achilles_df.quantile(0.1)
ceres_stats_df['q25_ceres'] = achilles_df.quantile(0.25)
ceres_stats_df['q50_ceres'] = achilles_df.quantile(0.5)
ceres_stats_df['q75_ceres'] = achilles_df.quantile(0.75)
ceres_stats_df['q90_ceres'] = achilles_df.quantile(0.9)
ceres_stats_df['min'] = achilles_df.min()
ceres_stats_df['max'] = achilles_df.max()
ceres_stats_df = ceres_stats_df.loc[ceres_stats_df.index.intersection(intersect_genes)]
#Save the dataframe as follows...
ceres_stats_df.to_parquet(f'{processed_path}/ceres_stats_df.parquet')
