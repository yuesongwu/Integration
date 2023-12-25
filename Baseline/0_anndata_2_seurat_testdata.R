library(Seurat)
library(SeuratDisk)
library(Matrix)
library(dplyr)
library(ggplot2)
# library(rhdf5)
library(reticulate)

setwd('/mnt/home/wuyueson/WorkFiles/scInformer/Baseline')

reticulate::use_condaenv("scInformer", required = TRUE)
scanpy =  import("scanpy")
hdf5plugin = import("hdf5plugin")

#TODO: to test on other dataset, change DATA_PATH
DATA_PATH = '../data/HLCA_zstd.h5ad'
data = scanpy$read(DATA_PATH)
##  create test dataset 
data = data[data$obs$study %in% c('Teichmann_Meyer_2019', 'Seibold_2020')]
count =  t(data$X) 
rownames(count) = data$var_names$to_list()
colnames(count) = data$obs_names$to_list()

merged_seurat_test = CreateSeuratObject(count)
merged_seurat_test = AddMetaData(merged_seurat_test, data$obs)
merged_seurat_test
saveRDS(merged_seurat_test, file = '../data/HLCA_test.rds' )




## 1) create full dataset 
# count =  t(data$X) 
# rownames(count) = data$var_names$to_list()
# colnames(count) = data$obs_names$to_list()

# merged_seurat = CreateSeuratObject(count)
# merged_seurat = AddMetaData(merged_seurat, data$obs)
# merged_seurat 
# SaveH5Seurat(merged_seurat, filename = '../data/HLCA.h5Seurat' )

