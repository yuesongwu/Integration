library(Seurat)
# library(SeuratDisk)
library(SeuratWrappers)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reticulate)
reticulate::use_condaenv("leiden", required = TRUE)
ad <- import("anndata")
scanpy <-  import("scanpy")
hdf5plugin <- import("hdf5plugin")
setwd('/mnt/home/wuyueson/WorkFiles/scInformer/Baseline')


args = commandArgs(T)
dataset = args[1]
batch = args[2]
test = as.logical(args[3])

print('FastMNN')
#---------- load anndata ----------------
if (test){
    merged_seurat = readRDS('../data/HLCA_test.rds')
    hvgs = read.csv('../data/hvg2000_test.csv',header = TRUE)$X
}else{
    DATA_PATH = paste0('../data/',dataset,'.h5ad')
    data = scanpy$read(DATA_PATH)

    count =  t(data$X) 
    rownames(count) = data$var_names$to_list()
    colnames(count) = data$obs_names$to_list()

    merged_seurat = CreateSeuratObject(count)
    merged_seurat = AddMetaData(merged_seurat, data$obs)
    hvgs = read.csv('../data/hvg2000.csv',header = TRUE)$X
}
print(merged_seurat)


## use the same 2000hvg for Harmony
## seurat v5
## tutorial: https://github.com/satijalab/seurat-wrappers/blob/master/R/fast_mnn_v5.R
merged_seurat[['RNA']] =  split(merged_seurat[["RNA"]], f = merged_seurat$batch)
merged_seurat <- merged_seurat %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst") %>% 
    ScaleData() %>%
    RunPCA()

## ------------- FastMNN ----------------
FastMNN_seurat <- IntegrateLayers(object = merged_seurat, method = FastMNNIntegration,
    new.reduction = 'mnn', verbose = FALSE,features = hvgs)

# saveRDS(FastMNN_seurat, paste0("integrated_data/",dataset,"_FastMNN.rds" ))
emb = Embeddings(FastMNN_seurat, 'mnn')
## save emb as anndata, for the downstream scib metrics 
embed_adata = ad$AnnData(X = emb, obs = data.frame(row.names = rownames(emb)) , var = data.frame(row.names = colnames(emb)))
embed_adata$obs = FastMNN_seurat@meta.data
embed_adata$write_h5ad(paste0("integrated_data/",dataset,"_embed_FastMNN.h5ad"))




