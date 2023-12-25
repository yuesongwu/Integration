# import scib
import anndata as ad
import torch
import scvi
import scanpy as sc
import hdf5plugin
import argparse



#------------------------------------------- set up scvi ------------------
torch.manual_seed(0)
scvi.settings.seed = 0
sc.settings.verbosity = 0
#------------------------------------------- set up parameters ------------------
parser = argparse.ArgumentParser(description="model", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--data',default="HLCA_zstd_2000",type = str)
parser.add_argument('--test',action='store_true', default=False)
parser.add_argument('--batch_correction',default="batch",type = str)
args = parser.parse_args()


Integrated_data_path = f"integrated_data/{args.data}_scvi.h5ad"
Embed_data_path = f"integrated_data/{args.data}_embed_scvi.h5ad"
model_save_path = f"trained_models/{args.data}_scvi"


#------------------------ load data ------------------
print('scvi')
combined_adata = sc.read('../data/'+args.data + '.h5ad') 
print(combined_adata)
scvi.model.SCVI.setup_anndata(
    combined_adata,
    layer="counts",
    categorical_covariate_keys=[args.batch_correction]
)


#! model train 
model = scvi.model.SCVI(combined_adata)
# model
model.train()
# model save and loading
model.save(model_save_path)
# model = scvi.model.SCVI.load("trained_models/scvi", adata=combined_adata, use_gpu=True)

## get latent embedding 
latent = model.get_latent_representation()
combined_adata.obsm["X_latent"] = latent
sc.pp.neighbors(combined_adata, use_rep="X_latent")
sc.tl.umap(combined_adata, min_dist=0.3)

## save integrated data
combined_adata.write_h5ad(Integrated_data_path,compression=hdf5plugin.FILTERS["zstd"])
## save embed 
embed_adata = ad.AnnData(X = combined_adata.obsm["X_latent"])
embed_adata.obs = combined_adata.obs
embed_adata
embed_adata.write_h5ad(Embed_data_path,compression=hdf5plugin.FILTERS["zstd"])

# sc.pl.umap(
#     combined_adata,
#     color=["cell_type"],palette="tab20",
#     frameon=False,save = (args.data+ "_scvi_latent_celltype.pdf")
# )
# sc.pl.umap(
#     combined_adata,
#     color=["batch"],palette="tab20",
#     ncols=2,
#     frameon=False,save =( args.data + "_scvi_latent_batch.pdf")
# )

