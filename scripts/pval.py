def gene_confidence(gene_dict):
    gene_confidence_dict = {}
    global model_genes
    confidence_zero_genes = set(model_genes) - set(gene_dict.keys())
    
    for gene in confidence_zero_genes:
        gene_confidence_dict[gene] = 0
        
    for k,v in gene_dict.items():
        gene_expr_stats_df = expression_stats_df.loc[k]
        q_1 = gene_expr_stats_df['q10_expr']
        q_2 = gene_expr_stats_df['q50_expr']
        q_3 = gene_expr_stats_df['q75_expr']
        q_4 = gene_expr_stats_df['q90_expr']
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