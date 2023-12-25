import warnings
import anndata as ad
import Cell_BLAST as cb
import scanpy as sc
import os
import hdf5plugin
import sys
import argparse
import pandas as pd
#------------------------------------------- set up cb ------------------
warnings.filterwarnings("ignore")
cb.config.N_JOBS = 1
cb.config.RANDOM_SEED = 0
#------------------------------------------- set up parameters ------------------
parser = argparse.ArgumentParser(description="model", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--data',default="HLCA_zstd",type = str)
parser.add_argument('--test',action='store_true', default=False)
parser.add_argument('--batch_correction',default="batch",type = str)
parser.add_argument('--hvgfrac', default=0.5, type=float)
parser.add_argument('--latentdim', default=10, type=int)
parser.add_argument('--catdim', default=None, type=int)
parser.add_argument('--seed', default=0, type=int)
args = parser.parse_args()

Integrated_data_path = f"integrated_data/{args.data}_cb.h5ad"
Embed_data_path = f"integrated_data/{args.data}_embed_cb.h5ad"
model_save_path = f"trained_models/{args.data}_cb"


#------------------------ load data and integrate ------------------
print('CellBlast')
combined_adata = sc.read('../data/'+ args.data + '.h5ad') 
print(combined_adata)
combined_adata.X = combined_adata.layers['counts']

# if args.test:
#     combined_adata = combined_adata[combined_adata.obs['study'].isin(['Teichmann_Meyer_2019', 'Seibold_2020'])].copy()
#     print(combined_adata)
#     Integrated_data_path = f"integrated_data/{args.data}_cb_test.h5ad"
#     model_save_path = 'trained_models/cb_test'
    
# axes = cb.data.find_variable_genes(combined_adata,  min_group_frac=args.hvgfrac)
# print(combined_adata.var["variable_genes"].sum()) #0.2-1172

## use same hvg as scvi
print("start training") 
# hvg_df = pd.read_csv('../data/hvg2000.csv')
model_tosave =cb.directi.fit_DIRECTi(
        combined_adata,
        batch_effect = args.batch_correction,
        latent_dim = args.latentdim, cat_dim=args.catdim, random_seed=args.seed
    )
# model_tosave =cb.directi.fit_DIRECTi(
#         combined_adata, genes = combined_adata.var.query("variable_genes").index,
#         batch_effect = args.batch_correction,
#         latent_dim = args.latentdim, cat_dim=args.catdim, random_seed=args.seed
#     )
#-------- save model--------------
model_tosave.save(model_save_path )

#-------- save data and embed -------------
model = model_tosave
del model_tosave
combined_adata.obsm["X_latent"] = model.inference(combined_adata)
sc.pp.neighbors(combined_adata, use_rep="X_latent")
sc.tl.umap(combined_adata)
combined_adata.write_h5ad(Integrated_data_path,compression=hdf5plugin.FILTERS["zstd"])

embed_adata = ad.AnnData(X = combined_adata.obsm["X_latent"])
embed_adata.obs = combined_adata.obs
embed_adata
embed_adata.write_h5ad(Embed_data_path,compression=hdf5plugin.FILTERS["zstd"])


