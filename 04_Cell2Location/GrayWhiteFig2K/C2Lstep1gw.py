# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AN1792 Project                              -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Written by: Thomas Watson
# Summary: Gray/White Cell2Location Step 1
#
#-------------------------------------------------------------------------------

import os
import scanpy as sc
import numpy as np
import pandas as pd
import scipy as sp
import torch
import cell2location as c2l
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib as mpl
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

# Due to class imbalance, we tried downsampling cell types in various reference datasets with more cells than this "target_cells" threshold.
# After switching to the ROSMAP reference, this downsampling is no longer necessary.
# This number is now set much larger than the number of cells in the entire reference so the "downsampling" steps do nothing.
target_cells = 500000

print("INFO: Started!")
all_sample_sc = sc.read_h5ad("C2L/data.h5ad")
gw_ref = sc.read_h5ad('C2L/gwreferencedata.h5ad')
results_folder = "C2L/results"

# remove mitochondria-encoded (MT) genes
all_sample_sc.var['MT_gene'] = [gene.startswith('MT-') for gene in all_sample_sc.var_names]

# remove MT genes for spatial mapping (keeping their counts in the object)
all_sample_sc.obsm['MT'] = all_sample_sc[:, all_sample_sc.var['MT_gene'].values].X.toarray()
all_sample_sc = all_sample_sc[:, ~all_sample_sc.var['MT_gene'].values]

# remove hemoglobin-encoded (HB) genes
all_sample_sc.var['HB_gene'] = [gene.startswith('HB') for gene in all_sample_sc.var_names]

# remove HB genes for spatial mapping (keeping their counts in the object)
all_sample_sc.obsm['HB'] = all_sample_sc[:, all_sample_sc.var['HB_gene'].values].X.toarray()
all_sample_sc = all_sample_sc[:, ~all_sample_sc.var['HB_gene'].values]

# This is gray/white - remove mng
all_sample_sc = all_sample_sc[all_sample_sc.obs["manual_annotation"] != "meninges"]

# Make sure data type is int for C2L
# https://github.com/vitkl/cell2location_paper/blob/master/notebooks/lymph_nodes_analysis/analysis_1_estimating_reference_cell_type_signatures.ipynb
all_sample_sc.X = all_sample_sc.X.astype('int')
gw_ref.X = gw_ref.X.astype('int')

# This is the downsampling bit mentioned earlier - target_cells = 500000 so it does nothing
gw_ref_split = [gw_ref[gw_ref.obs["ct"]==clust] for clust in gw_ref.obs["ct"].unique()]
for dat in gw_ref_split:
    if dat.n_obs > target_cells:
        sc.pp.subsample(dat, n_obs=target_cells)

# put data back together
gw_ref_downsampled = gw_ref_split[0].concatenate(*gw_ref_split[1:])
gw_ref_downsampled.var['SYMBOL'] = gw_ref_downsampled.var.index

# Ensure same gene sets
shared_features = [
    feature for feature in all_sample_sc.var_names if feature in gw_ref_downsampled.var_names
]

# Make sure plot device is off
plt.close()

# Subset to shared genes
all_sample_sc = all_sample_sc[:, shared_features].copy()
gw_ref_downsampled = gw_ref_downsampled[:, shared_features].copy()

# permissive gene filtering per C2L's recommendation
selected = filter_genes(gw_ref_downsampled, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# save plot
plt.savefig(f"{results_folder}/plots/filter_genes_1223.pdf", bbox_inches = "tight")
plt.close()

# filter the object
gw_ref_downsampled = gw_ref_downsampled[:, selected].copy()

# set up NB regression
c2l.models.RegressionModel.setup_anndata(adata=gw_ref_downsampled,
                                         batch_key='batch',
                                         labels_key='ct'
                                         )
                                         
# NB reg
mod = RegressionModel(gw_ref_downsampled)
# view anndata_setup as a sanity check
mod.view_anndata_setup()
print("INFO: Training!")
# train model
mod.train(max_epochs=500, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)

# Plot training loss - plotting device has issues with the HPC this is run on
# We confirmed convergence of loss in the log files
mod.plot_history(500)
plt.savefig(f"{results_folder}/plots/training_loss_1223.pdf")
plt.close()

# Export posterior
gw_ref_downsampled = mod.export_posterior(
    gw_ref_downsampled, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in gw_ref_downsampled.varm.keys():
    inf_aver = gw_ref_downsampled.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' 
                                    for i in gw_ref_downsampled.uns['mod']['factor_names']]].copy()
else:
    inf_aver = gw_refdownsampled.var[[f'means_per_cluster_mu_fg_{i}' 
                                    for i in gw_ref_downsampled.uns['mod']['factor_names']]].copy()
inf_aver.columns = gw_ref_downsampled.uns['mod']['factor_names']

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(gw_ref_downsampled.var_names, inf_aver.index)
all_sample_sc = all_sample_sc[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# save predictions to CSV
inf_aver.to_csv(f"{results_folder}/inf_aver_gw_1223.csv")

# save plot - the plot device on HPC prints both plots in the same space on top of each other!! 
# github issue is open but unresolved
# https://github.com/BayraktarLab/cell2location/issues/341
mod.plot_QC()
plt.savefig(f"{results_folder}/plots/training_qc_1223.pdf")
plt.close()

# print to log file
print("INFO: Saving!")

# Save model
mod.save(f"{results_folder}/models/model_1223", overwrite=True)

# Save anndata objects with results
adata_file = f"{results_folder}/adata_objects/gw_ref_1223.h5ad"
gw_ref_downsampled.write(adata_file)
sc_file = f"{results_folder}/adata_objects/all_sample_sc_after_training_1223.h5ad"
all_sample_sc.write(sc_file)
