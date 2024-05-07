# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AN1792 Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Written by: Thomas Watson
# Summary: Save C2L predicted abundances
#
#-------------------------------------------------------------------------------

import os
import warnings
warnings.filterwarnings('ignore')
import torch
import pickle
import sklearn
import scipy as sp
import numpy as np
import scanpy as sc
import pandas as pd
import squidpy as sq
import anndata as ad
import matplotlib as mpl
import cell2location as c2l
import matplotlib.pyplot as plt
from matplotlib import rcParams
from cell2location.plt import plot_spatial
from cell2location import run_colocation
from cell2location.utils import select_slide
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from sklearn.model_selection import train_test_split
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

samples = ["A34285", "A34583.2", "A34644", "A34717", "A34933.2", "A34995", "A35038", "AN1792.102.10",
           "AN1792.102.11", "AN1792.102.15", "AN1792.102.16", "AN1792.102.17", "AN1792.102.19", "AN1792.102.20", "AN1792.102.21", 
           "AN1792.102.22", "AN1792.102.4",  "AN1792.102.6",  "AN1792.102.7",  "AN1792.102.8",  "AN1792.102.1",  "N35127N",
           "A18.148", "A34291A", "A34992A", "AX21.92"]

for sample in samples:
    
    adatas = []
    spatial_dictm = {}
    spatial_dictg = {}
    spatial_dictc = {}

    gw = sc.read_h5ad(f"/{sample}_1120.h5ad")
    mng = sc.read_h5ad(f"/{sample}_120.h5ad")
    dat2 = sc.read_visium(f"{sample}/outs/", 
                      count_file='filtered_feature_bc_matrix.h5', 
                      source_image_path = 'spatial/tissue_lowres_image.png')

    dat2.obs['barcode'] = dat2.obs.index
    
    model_folder_gw = f"/models/{sample}/"
    model_gw = c2l.models.Cell2location.load(model_folder_gw, gw)

    expected_dict_gw = model_gw.module.model.compute_expected_per_cell_type(
        gw.uns['mod']["post_sample_q05"], model_gw.adata_manager
    )

    for i, n in enumerate(model_gw.factor_names_):
        gw.layers[n] = expected_dict_gw['mu'][i]
    
    sc.pp.neighbors(gw, use_rep="q05_cell_abundance_w_sf")
    sc.tl.leiden(gw, resolution=0.5)
    gw.obs["region_cluster"] = gw.obs["leiden"].astype("category")
    sc.tl.umap(gw, min_dist=0.3, spread=1)
    
    model_folder_mng = f"/models/{sample}/"
    model_mng = c2l.models.Cell2location.load(model_folder_mng, mng)


    expected_dict_mng = model_mng.module.model.compute_expected_per_cell_type(
      mng.uns['mod']["post_sample_q05"], model_mng.adata_manager
    )

    for i, n in enumerate(model_mng.factor_names_):
        mng.layers[n] = expected_dict_mng['mu'][i]
        
    sc.pp.neighbors(mng, use_rep="q05_cell_abundance_w_sf")
    sc.tl.leiden(mng, resolution=0.5)
    mng.obs["region_cluster"] = mng.obs["leiden"].astype("category")
    sc.tl.umap(mng, min_dist=0.3, spread=1)

    dat2m = dat2[dat2.obs['barcode'].isin(mng.obs['barcode'])]

    spatial_dictm.update(dat2m.uns["spatial"])

    mng.uns['spatial'] = spatial_dictm
    mng.obsm['spatial'] = dat2m.obsm['spatial'].astype("int")

    dat2g = dat2[dat2.obs['barcode'].isin(gw.obs['barcode'])]

    spatial_dictg.update(dat2g.uns["spatial"])

    gw.uns['spatial'] = spatial_dictg
    gw.obsm['spatial'] = dat2g.obsm['spatial'].astype("int")

    merged = ad.concat([gw, mng], join = "outer")
    
    adatas.append(merged)

    combined = ad.concat(adatas, join = "outer")

    dat2c = dat2[dat2.obs['barcode'].isin(combined.obs['barcode'])]

    spatial_dictc.update(dat2c.uns["spatial"])

    combined.uns['spatial'] = spatial_dictc
    
    combined.obs.to_csv(f"/{sample}_combined_meta.csv")
