import pandas as pd

ceres_stats_df = pd.read_parquet('../data/processed/CRISPR_gene_effect_processed.parquet')
sbml_fname = './models/Recon3DModel_301_no_isoforms_gene_symbols_rpmi_trimed.xml'

model = read_sbml_model(sbml_fname)
