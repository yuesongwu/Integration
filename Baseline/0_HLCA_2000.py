# import scib
import anndata as ad
import scanpy as sc
import hdf5plugin
import argparse
import pandas as pd

#------------------------------------------- set up parameters ------------------
parser = argparse.ArgumentParser(description="model", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--data',default="HLCA_zstd",type = str)
parser.add_argument('--test',action='store_true', default=False)
parser.add_argument('--batch_correction',default="batch",type = str)
args = parser.parse_args()

print('scvi')
combined_adata = sc.read('../data/'+args.data + '.h5ad') 
print(combined_adata)
if args.test:
    combined_adata = combined_adata[combined_adata.obs['study'].isin(['Teichmann_Meyer_2019', 'Seibold_2020'])].copy()
    print(combined_adata)
    Integrated_data_path = f"integrated_data/{args.data}_scvi_test.h5ad"
    model_save_path = 'trained_models/scvi_test'
# combined_adata.obs.cell_type
# combined_adata.obs.batch
sc.pp.filter_genes(combined_adata, min_counts=10)
# print(combined_adata.X[0:10,0:10])

combined_adata.layers["counts"] = combined_adata.X.copy()  # preserve counts
sc.pp.normalize_total(combined_adata, target_sum=1e4)
sc.pp.log1p(combined_adata)
combined_adata.raw = combined_adata # freeze the state in `.raw`


sc.pp.highly_variable_genes(
    combined_adata,
    n_top_genes=2000,
    subset=True,
    layer="counts",
    flavor="seurat_v3"
)
print(combined_adata)
combined_adata.write_h5ad(f"../data/{args.data}_2000.h5ad",compression=hdf5plugin.FILTERS["zstd"])
pd.DataFrame({'hvg':combined_adata.var['highly_variable']}).to_csv('../data/hvg2000.csv')