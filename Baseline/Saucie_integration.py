import sklearn.decomposition
import scanpy as sc
import hdf5plugin
from scipy.sparse import csr_matrix
import anndata as ad
import sys
import SAUCIE
import os
import argparse
import h5py

#------------------------------------------- set up parameters ------------------
parser = argparse.ArgumentParser(description="model", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--data',default="HLCA_zstd_2000",type = str)
parser.add_argument('--test',action='store_true', default=False)
parser.add_argument('--batch_correction',default="batch",type = str)
args = parser.parse_args()

Integrated_data_path = f"integrated_data/{args.data}_saucie.h5ad"
Embed_data_path = f"integrated_data/{args.data}_embed_saucie.h5ad"
model_save_path = f"trained_models/{args.data}_saucie"

#------------------------------------------- load by h5py manually------------------
# file = h5py.File('/mnt/home/wuyueson/WorkFiles/scInformer/data/HLCA_zstd_2000.h5ad')
# # # file = h5py.File('/mnt/gs21/scratch/wuyueson/CellBlast/cellblast_input/HLCA_nozstd.h5ad')
# # # file = h5py.File('/mnt/gs21/scratch/wuyueson/CellBlast/cellblast_input/benchmark_inner_nozstd.h5ad')
# list(file.keys())
# list(file['layers']['counts'].keys())
# # list(file['X'].keys())
# # list(file['obs'].keys())
# # list(file['obs']['cell_type'].keys())

# # #csr_matrix((data, indices, indptr), [shape=(M, N)])
# count = csr_matrix((file['layers']['counts']['data'],file['layers']['counts']['indices'], file['layers']['counts']['indptr'] ))
# print(count[0:10,0:10])
# adata = anndata.AnnData(X = count)
# adata.var_names =  [ str(x)[2:-1] for x in  file['var']['_index']] #HLCA
# adata.obs_names = [ str(x)[2:-1] for x in  file['obs']['_index']]


# def get_labels(labelname):
#     #list() to retrieve data in the 'Dataset' object of  <file['obs']['cell_type']['categories']>
#     temp = [ str(x)[2:-1] for x in  file['obs'][labelname]['categories']]
#     ind_temp = list(file['obs'][labelname]['codes'])
#     out = [temp[i] for i in ind_temp]
#     return out

# adata.obs['cell_type'] = get_labels('cell_type')
# adata.obs['study'] = get_labels('study')
# adata.obs['batch'] = get_labels('batch')

# # adata.obs['cell_type'].value_counts()
# # adata.obs['study'].value_counts()
# # adata.obs['batch'].value_counts()
# adata 
# adata.write_h5ad( f"../data/{args.data}_saucie.h5ad",compression=hdf5plugin.FILTERS["zstd"])



#------------------------------
print('saucie')
## tutorial: https://github.com/KrishnaswamyLab/SAUCIE/tree/master?tab=readme-ov-file    example
# ## https://colab.research.google.com/github/KrishnaswamyLab/SingleCellWorkshop/blob/master/exercises/Deep_Learning/notebooks/02_Answers_Exploratory_analysis_of_single_cell_data_with_SAUCIE.ipynb
adata = sc.read(  f"../data/{args.data}_saucie.h5ad")
batch = args.batch_correction

pca_op = sklearn.decomposition.PCA(100)
if isinstance(adata.X, csr_matrix):
    expr = adata.X.A
else:
    expr = adata.X
data = pca_op.fit_transform(expr)
print("saucie start")
#lambda_b for batch correction
saucie = SAUCIE.SAUCIE(100, lambda_b=0.5)
loader_train = SAUCIE.Loader(data, labels=adata.obs[batch].cat.codes, shuffle=True)
loader_eval = SAUCIE.Loader(data, labels=adata.obs[batch].cat.codes, shuffle=False)
print("saucie train")
saucie.train(loader_train, steps=5000)
ret = adata.copy()
ret.obsm["X_latent"] = saucie.get_reconstruction(loader_eval)[0]
ret.X = pca_op.inverse_transform(ret.obsm["X_latent"])

# pca umap 
sc.tl.pca(ret)
sc.pp.neighbors(ret, n_pcs=30, n_neighbors=20)
sc.tl.umap(ret)


## save integrated data
ret.write_h5ad(Integrated_data_path,compression=hdf5plugin.FILTERS["zstd"])
## save embed 
embed_adata = ad.AnnData(X = ret.obsm["X_latent"])
embed_adata.obs = ret.obs
print(embed_adata)
embed_adata.write_h5ad(Embed_data_path,compression=hdf5plugin.FILTERS["zstd"])

