import pandas as pd
import re

raw_path = './data/raw'
processed_path = './data/processed'

df_sample_info = pd.read_csv(f'{raw_path}/sample_info.csv', index_col=0)
ccle_names = df_sample_info.to_dict()['CCLE_Name']

df_crispr_gene = pd.read_csv(f'{raw_path}/CRISPR_gene_effect.csv', index_col=0)

df_crispr_gene= df_crispr_gene.rename(ccle_names, axis=0)
column_rename = {elem: re.sub(r' \([0-9]+\)$',r'', elem) for elem in df_crispr_gene.columns}
df_crispr_gene = df_crispr_gene.rename(column_rename, axis=1)
df_crispr_gene = df_crispr_gene.fillna(0)
print(df_crispr_gene.head())

df_crispr_gene.to_parquet(f'{processed_path}/CRISPR_gene_effect_processed.parquet')