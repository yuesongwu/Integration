import scanpy as sc
import pandas as pd
import hdf5plugin
import os
import scanorama
import anndata
import argparse
import anndata as ad



#------------------------------------------- set up parameters ------------------
parser = argparse.ArgumentParser(description="model", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--data',default="HLCA_zstd",type = str)
parser.add_argument('--test',action='store_true', default=False)
parser.add_argument('--batch_correction',default="batch",type = str)
args = parser.parse_args()

Integrated_data_path = f"integrated_data/{args.data}_scanorama.h5ad"
Embed_data_path = f"integrated_data/{args.data}_embed_scanorama.h5ad"
model_save_path = f"trained_models/{args.data}_scanorama"


# load 
print('scanorama')
combined_adata = sc.read('../data/'+args.data + '.h5ad') 
print(combined_adata)

# split per batch into new objects
split = []
batch_categories = combined_adata.obs[args.batch_correction].cat.categories

for i in batch_categories:
    split.append(combined_adata[combined_adata.obs[args.batch_correction] == i].copy())

corrected = scanorama.correct_scanpy(split,  return_dimred=True)
merged = anndata.AnnData.concatenate(*corrected ,batch_key=args.batch_correction, batch_categories= batch_categories)

merged.obsm["X_latent"] = merged.obsm["X_scanorama"]
sc.pp.neighbors(merged, use_rep = "X_latent")
sc.tl.umap(merged)
merged.write_h5ad(Integrated_data_path,compression=hdf5plugin.FILTERS["zstd"])
embed_adata = ad.AnnData(X = merged.obsm["X_latent"])
embed_adata.obs = merged.obs
embed_adata
embed_adata.write_h5ad(Embed_data_path,compression=hdf5plugin.FILTERS["zstd"])








