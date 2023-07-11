import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from tqdm import tqdm 
import os


section_ids = ['151673','151674','151675','151676']


for section_id in tqdm(section_ids[0]):
    print(section_id)
    input_dir = "../data/DLPFC/data/{}/".format(section_id)
    adata = sc.read_visium(path=input_dir, count_file=section_id + '_filtered_feature_bc_matrix.h5', load_images=True)
    adata.var_names_make_unique(join="++")

    # read the annotation
    Ann_df = pd.read_csv(os.path.join(input_dir, section_id + '_truth.txt'), sep='\t', header=None, index_col=0)
    Ann_df.columns = ['Ground Truth']
    Ann_df[Ann_df.isna()] = "unknown"
    adata.obs['Ground Truth'] = Ann_df.loc[adata.obs_names, 'Ground Truth'].astype('category')
    
    # make spot name unique
    adata.obs_names = [x+'_'+section_id for x in adata.obs_names]
    
    # Constructing the spatial network
    # STAligner.Cal_Spatial_Net(adata, rad_cutoff=150) # the spatial network are saved in adata.uns[‘adj’]
    # # STAligner.Stats_Spatial_Net(adata) # plot the number of spatial neighbors
    
    # # Normalization
    # sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=5000)
    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)
    # adata = adata[:, adata.var['highly_variable']]

    # adj_list.append(adata.uns['adj'])
    
