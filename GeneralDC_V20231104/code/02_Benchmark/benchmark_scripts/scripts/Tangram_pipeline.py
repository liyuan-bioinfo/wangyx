import pandas as pd
import scanpy as sc
import numpy as np
import tangram as tg
import sys

sc_file_path = sys.argv[1]
spatial_file_path = sys.argv[2]
celltype_key = sys.argv[3]
output_file_path = sys.argv[4]


sc_adata = sc.read_h5ad(sc_file_path)
sp_adata = sc.read_h5ad(spatial_file_path)


# run Tangram
intersect = np.intersect1d(sc_adata.var_names, sp_adata.var_names)        
sp_adata = sp_adata[:, intersect].copy()
sc_adata = sc_adata[:, intersect].copy()
## Find DEG for sc
sc.pp.log1p(sc_adata)
sc.tl.rank_genes_groups(sc_adata, groupby=celltype_key, use_raw=False)

markers_df = pd.DataFrame(sc_adata.uns["rank_genes_groups"]["names"]).iloc[0:2000, :]

genes_sc = np.unique(markers_df.melt().value.values)
genes_st = sp_adata.var_names.values
genes = list(set(genes_sc).intersection(set(genes_st)))
tg.pp_adatas(sc_adata, sp_adata, genes=genes)

ad_map = tg.map_cells_to_space(
                sc_adata,
                sp_adata,
                mode='clusters',
                cluster_label=celltype_key)

tg.project_cell_annotations(ad_map, sp_adata, annotation=celltype_key)

celltype_density = sp_adata.obsm['tangram_ct_pred']
celltype_density = (celltype_density.T/celltype_density.sum(axis=1)).T

celltype_density.to_csv(f"{output_file_path}.csv")   

