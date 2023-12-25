# Integration


## Benchmark
for immune 
https://openproblems.bio/results/batch_integration_embed/

## Journal/Conference
bioinformatics:
https://academic.oup.com/bioinformatics/
https://academic.oup.com/bioinformatics/pages/submission_online
    https://academic.oup.com/bioinformatics/pages/instructions_for_authors
    https://academic.oup.com/bioinformatics/pages/submission_online
    
ISMB 2024:
https://www.iscb.org/ismb2024/home

## Baseline 
### Installation
git clone https://github.com/C0nc/scInformer.git
cd scInformer
conda create -n leiden python=3.9
conda activate leiden
conda install cudatoolkit=11.7 -c pytorch -c nvidia
conda install -c conda-forge faiss-gpu
pip install torch==1.13.1+cu117 torchvision==0.14.1+cu117 torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cu117
pip install --extra-index-url https://pypi.nvidia.com cudf-cu11==23.4.1 dask-cudf-cu11==23.4.1 cuml-cu11==23.4.1 cugraph-cu11==23.4.1 cucim==23.4.1
pip install einops ipdb pydance torchmetrics rapids-singlecell scvi-tools scib wandb hdf5plugin
pip install scib-metrics 
pip install captum
pip install scanorama 
pip install Cell-BLAST



conda create -n saucie python=3.6
conda activate saucie
git clone https://github.com/KrishnaswamyLab/SAUCIE.git
pip install -r SAUCIE/requirements.txt
pip install packaging
pip install natsort
pip install numba
pip install tqdm

R4.3.0; Seurat V5
FastMNN; Harmony

### Baseline HPCC Notes 
Baseline(GPU): scvi, Cell-BLAST   (on hpcc: request v100:gpu, b/c k80 is cuda11.4, too old)

scvi(raw) 0.5h; CB(raw) 1h; scanorama(lognorm) 1.5h; Harmony(lognorm) 1.5h(Seurat v4); FastMNN(lognorm) 13.5h (Seurat v5)
saucie 0.5h?

#TODO: need to rerun FastMNN, using same 2000 hvg (12/24 using 3000); try Harmony v5 (also 2000)

## Evaluation

