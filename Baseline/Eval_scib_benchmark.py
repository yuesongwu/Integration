import scanpy as sc
import hdf5plugin
import scib
import faiss
from scib_metrics.nearest_neighbors import NeighborsOutput
from scib_metrics.benchmark import Benchmarker,BioConservation
import argparse
import os
os.chdir('/mnt/home/wuyueson/WorkFiles/scInformer/Baseline')
import pandas as pd
import numpy as np
import pickle
# import rapids_singlecell as rsc


#--------------- functions -----------------

def faiss_hnsw_nn(X: np.ndarray, k: int):
    """Gpu HNSW nearest neighbor search using faiss.

    See https://github.com/nmslib/hnswlib/blob/master/ALGO_PARAMS.md
    for index param details.
    """
    X = np.ascontiguousarray(X, dtype=np.float32)
    res = faiss.StandardGpuResources()
    M = 32
    index = faiss.IndexHNSWFlat(X.shape[1], M, faiss.METRIC_L2)
    gpu_index = faiss.index_cpu_to_gpu(res, 0, index)
    gpu_index.add(X)
    distances, indices = gpu_index.search(X, k)
    del index
    del gpu_index
    # distances are squared
    return NeighborsOutput(indices=indices, distances=np.sqrt(distances))


def faiss_brute_force_nn(X: np.ndarray, k: int):
    """Gpu brute force nearest neighbor search using faiss."""
    X = np.ascontiguousarray(X, dtype=np.float32)
    res = faiss.StandardGpuResources()
    index = faiss.IndexFlatL2(X.shape[1])
    gpu_index = faiss.index_cpu_to_gpu(res, 0, index)
    gpu_index.add(X)
    distances, indices = gpu_index.search(X, k)
    del index
    del gpu_index
    # distances are squared
    return NeighborsOutput(indices=indices, distances=np.sqrt(distances))

#------------------------------------------- set up parameters ------------------
parser = argparse.ArgumentParser(description="model", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--data',default="HLCA_zstd_2000",type = str)
parser.add_argument('--batch_correction',default="batch",type = str)
args = parser.parse_args()
#tutorial: https://scib-metrics.readthedocs.io/en/stable/notebooks/large_scale.html

#------------------------ load data ------------------
adata = sc.read('../data/'+ args.data + '.h5ad')
sc.tl.pca(adata, n_comps=30, use_highly_variable=True)
adata.obsm['Unintegrated'] = adata.obsm["X_pca"]

## out: pd
Embed_path = './integrated_data/'
files = [x for x in os.listdir(Embed_path) if '_embed_' in x ]
# eval_df = pd.DataFrame(index = ['nmi','ari','cell_cycle'],columns =[f.split('_embed_')[1].replace('.h5ad','') for f in files] )
res_list = []
method_list = ["Unintegrated"]
for file in files[0:1]:
    adata_embedding = sc.read(Embed_path + file)
    method = file.split('_embed_')[1].replace('.h5ad','')
    print(method)
    adata.obsm[method] = adata_embedding.X
    method_list.append(method)
    
biocons = BioConservation(isolated_labels=False)
bm = Benchmarker(
    adata,
    batch_key="batch",
    label_key="cell_type",
    embedding_obsm_keys=["Unintegrated", "scanorama"],
    pre_integrated_embedding_obsm_key="X_pca",
    bio_conservation_metrics=biocons,
    n_jobs=-1,
)
bm.prepare(neighbor_computer=faiss_brute_force_nn) #! cudaMalloc error out of memory [2]
bm.benchmark()




    
with open('evaluation/metric_ari_test.pkl','rb') as f:
    pickle.dump(res_list,f)

# scib.me.metrics_fast()
# sc.tl.pca(adata, n_comps=30, use_highly_variable=True)
# adata.obsm["Unintegrated"] = adata.obsm["X_pca"]
# print(adata)

# bm = Benchmarker(
#     adata,
#     batch_key="batch",
#     label_key="cell_type",
#     embedding_obsm_keys= method_list,
#     n_jobs=6,
# )
# bm.benchmark()

# with open('evaluation/metric_test.pkl','rb') as f:
#     pickle.dump(bm,f)
# bm.plot_results_table()




    