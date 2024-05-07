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
# Summary: Gray/White Cell2Location Step 2
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

print("INFO: Started!")
results_folder = "C2L/results"
run_name = f'{results_folder}/cell2location_map'

dat = sc.read_h5ad("gw/all_sample_sc_after_training_1223.h5ad")
gw_ref_downsampled = sc.read_h5ad("gw_ref_1223.h5ad")

for sample in dat.obs["sample_id"].unique():
  
  all_sample_sc = dat[dat.obs["sample_id"] == sample]
  
  all_sample_sc = all_sample_sc[all_sample_sc.obs["manual_annotation"] != "meninges"]

  # load training posterior probabilities
  inf_aver = pd.read_csv("/C2Lresults/inf_aver_gw.csv")

  # index vals saved in column 1 but not index so make sure inf_aver.index = gw_ref_downsamled.var_names before running the model
  inf_aver.index = inf_aver.iloc[:]['Unnamed: 0'].astype("str")
  inf_aver.drop(columns=inf_aver.columns[0], axis=1,  inplace=True)
  all_sample_sc.var_names = all_sample_sc.var_names.astype("str")
  
  # prepare anndata for cell2location model
  c2l.models.Cell2location.setup_anndata(adata=all_sample_sc, batch_key= "sample_id")

  mod = c2l.models.Cell2location(
      all_sample_sc, cell_state_df=inf_aver,
      # the expected average cell abundance: tissue-dependent 
      # hyper-prior which can be estimated from paired histology:
      N_cells_per_location=7,
      # hyperparameter controlling normalisation of
      # within-experiment variation in RNA detection:
      detection_alpha=20
  ) 

  mod.train(max_epochs=5000, 
            # train using full data (batch_size=None)
            batch_size=None, 
            # use all data points in training because 
            # we need to estimate cell abundance at all locations
            train_size=1,
            use_gpu=True,
           )

  # plot ELBO loss history during training, removing first 100 epochs from the plot
  plt.close()
  mod.plot_history(5000)
  plt.legend(labels=['full data training'])
  plt.savefig(f"{results_folder}plots/{sample}_training_loss_1223.pdf", bbox_inches = "tight")
  plt.close()
  ###########################################################################################################
  all_sample_sc = mod.export_posterior(
      all_sample_sc, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
  )

  # Save model
  mod.save(f"{results_folder}models/{sample}_1223", overwrite=True)
  plt.close()
  mod.plot_QC()
  plt.savefig(f"{results_folder}plots/{sample}_training_QC_1223.pdf", bbox_inches = "tight")
  plt.close()
  #############################################################################################################

  all_sample_sc.obs[all_sample_sc.uns['mod']['factor_names']] = all_sample_sc.obsm['q05_cell_abundance_w_sf']
  
  # expected_dict = mod.module.model.compute_expected_per_cell_type(mod.samples["post_sample_q05"], mod.adata_manager)
  # 
  # for i, n in enumerate(mod.factor_names_):
  #   all_sample_sc.layers[n] = expected_dict["mu"][i]
  
  ###########################################################################################################
  # Save anndata object with results
  adata_file = f"{results_folder}adata_objects/{sample}_1223.h5ad"
  
  all_sample_sc.write(adata_file)
  
  print(f"{sample} data saved")
  
else: 
    print("All data saved")
