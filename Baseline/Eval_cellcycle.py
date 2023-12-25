import scanpy as sc
import hdf5plugin
import scib
# from scib_metrics.benchmark import Benchmarker,BioConservation
import argparse
import os
os.chdir('/mnt/home/wuyueson/WorkFiles/scInformer/Baseline')
import pandas as pd
import pickle
# import rapids_singlecell as rsc
# import numpy as np

#------------------------------------------- set up parameters ------------------
parser = argparse.ArgumentParser(description="model", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--data',default="HLCA_zstd_2000",type = str)
parser.add_argument('--batch_correction',default="batch",type = str)
args = parser.parse_args()

##TODO: tutorial
#https://scib-metrics.readthedocs.io/en/stable/notebooks/large_scale.html

# bm = Benchmarker(
#     adata,
#     batch_key="batch",
#     label_key="cell_type",
#     embedding_obsm_keys=["Unintegrated", "Scanorama", "LIGER", "Harmony", "scVI", "scANVI"],
#     n_jobs=6,
# )
# bm.benchmark()

#------------------------ load data ------------------
adata = sc.read('../data/'+args.data + '.h5ad')
adata_unintegrated = adata.copy()
# sc.tl.pca(adata, n_comps=30, use_highly_variable=True)
# adata.obsm["Unintegrated"] = adata.obsm["X_pca"]



## out: pd
Embed_path = './integrated_data/'
files = [x for x in os.listdir(Embed_path) if '_embed_' in x ]
eval_df = pd.DataFrame(index = ['nmi','ari','cell_cycle'],columns =[f.split('_embed_')[1].replace('.h5ad','') for f in files] )
method_list = ["Unintegrated"]
res_list = []
for file in files:
    adata_embedding = sc.read(Embed_path + file)
    method = file.split('_embed_')[1].replace('.h5ad','')
    print(method)
    adata.obsm[method] = adata_embedding.X
    method_list.append(method)
    eval_df.at['cell_cycle',method] = scib.me.cell_cycle(adata_pre=adata_unintegrated, adata_post=adata,batch_key= args.batch_correction ,embed  = method,organism='human')
    print(eval_df) 
    res = scib.metrics.metrics(adata=adata_unintegrated, adata_int=adata,
                            batch_key  = args.batch_correction, label_key = 'cell_type',embed = method,
                            cell_cycle_ = True,organism='human') 
    # print(res)
    #ari_ = True,nmi_ = True, silhouette_ = True, cell_cycle_ = True,pcr_ = True
    # print(res)
    # res_list.append(res)
# print(eval_df)  

# with open('./evaluation/metric_cellcycle_test.pkl','wb') as f:
#     pickle.dump(eval_df,f)

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




    