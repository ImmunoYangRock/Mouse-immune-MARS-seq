###############################################################################
# Fig. S6d — human sepsis scRNA-seq (GSE279452): preprocessing (Python part)
#
# Covered panel:
#   Fig. S6d: PDCD4 expression in CD8 T cells across sepsis groups, and its
#             correlation with effector / naive markers (R part in
#             FigureS6d_human_sepsis.R)
#
# Workflow:
#   1. read 10x matrices of all samples (GSE279452_RAW)
#   2. QC outliers (genes 200-5000, MT% <= 10, HB% <= 0.1)
#   3. scanpy standard pipeline (HVG 3000, regress out MT%, PCA, leiden 0.8)
#   4. save processed object; convert to .rds for the R plotting part
#
# Source notebook: lps_human.ipynb
###############################################################################

import os
import re
import numpy as np
import pandas as pd
import scanpy as sc  # 1.11.4
import anndata as ad  # 0.11.3
import seaborn as sns
from matplotlib import pyplot as plt

# os.chdir('lps_human')  # folder containing GSE279452_RAW/

# =============================================================================
# 1. read all samples
# =============================================================================
# !unzip GSE279452_RAW.zip
## each sample folder contains barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz
file_list = [i for i in os.listdir('./GSE279452_RAW/') if 'gz' not in i]

adata_list = []
for file in file_list:
    adata = sc.read_10x_mtx(f'./GSE279452_RAW/{file}/')
    adata.obs['sample'] = file.split('_')[1]
    adata.obs_names = [file.split('_')[1] + '_' + i for i in adata.obs_names]
    adata_list.append(adata)
adata = ad.concat(adata_list)

# =============================================================================
# 2. QC
# =============================================================================
def get_qc_outliers(adata):
    ## MT / Ribo / HB gene flags
    adata.var["MT"] = adata.var_names.str.startswith("MT-")
    adata.var["RB"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["HB"] = adata.var_names.str.contains(r"^HB[ABDEGMQZ]\d*(?!\w)")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["MT", "RB", "HB"], inplace=True,
                               log1p=True, percent_top=[20])
    ## outlier thresholds
    adata.obs["GE_outlier"] = (adata.obs["n_genes_by_counts"] < 200) | (adata.obs["n_genes_by_counts"] > 5000)
    adata.obs["MT_outlier"] = adata.obs["pct_counts_MT"] > 10
    adata.obs["HB_outlier"] = adata.obs["pct_counts_HB"] > 0.1
    print('GENE outliers: ', adata.obs.GE_outlier.value_counts())
    print('MT% outliers: ', adata.obs.MT_outlier.value_counts())
    print('HB% outliers: ', adata.obs.HB_outlier.value_counts())
    return adata


def remove_qc_outliers(adata):
    print(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.GE_outlier) & (~adata.obs.MT_outlier) & (~adata.obs.HB_outlier)].copy()
    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
    return adata


get_qc_outliers(adata)
adata = remove_qc_outliers(adata)

# =============================================================================
# 3. scanpy standard pipeline (leiden 0.8)
# =============================================================================
def scanpy_pip(adata, pcs=50):
    ## filter genes expressed in less than 5 cells
    sc.pp.filter_genes(adata, min_cells=5)
    ## normalize to median total counts, log1p
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    ## feature selection: top 3000 HVGs (batch-aware)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key='sample')
    ## remove MT / LINC / HB / IG / TCR from HVGs
    bad_features = adata.var_names[adata.var_names.str.contains(
        r'^(MT-|MTMNR|RP[SL]|LINC|^IG[HKL]|KRT\d+|^TR[AB][VJ]|^HB[ABDEGMQZ]\d*(?!\w))', regex=True)]
    for i in adata.var.highly_variable.index:
        if i in bad_features:
            adata.var.highly_variable[i] = False
    print('HVG resident numbers: ', adata.var.highly_variable.sum())
    ## freeze the log1p-normalized layer for visualization
    adata.raw = adata.copy()
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['pct_counts_MT'])
    sc.pp.scale(adata, max_value=10)
    ## PCA, neighbors, clustering, UMAP
    sc.pp.pca(adata, n_comps=50)
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
    sc.pp.neighbors(adata, n_pcs=pcs)
    sc.tl.leiden(adata, resolution=0.8, key_added='leiden_8')
    sc.tl.umap(adata)
    return adata


adata = scanpy_pip(adata, pcs=35)

## group = sample name without digits (e.g. HC, IC, Abd-S, Res-S, Skin-S, Uri-S,
## CNS-S, PHC, PIC, Abd-PS, Res-PS)
adata.obs['group'] = [re.sub(r'\d', '', i) for i in adata.obs['sample']]

## QC / annotation views
sc.pl.umap(adata, color=['sample', 'group', 'leiden_8'], cmap='RdPu', ncols=6)
sc.pl.umap(adata, color=['PTPRC', 'CD3D', 'CD19', 'CD4', 'CD8A', 'LYZ', 'ITGAM', 'C1QA',
                         'S100A8', 'S100A9', 'CD74', 'CCR6', 'RORC', 'GATA3', 'NCAM1',
                         'EOMES', 'PDCD4'], cmap='RdPu', ncols=6)

## lineage markers for cluster annotation (leiden_8 == '0' = CD8 T cells)
GD = {'IMMUNE': ['PTPRC'],
      'T': ['CD3D', 'CD3E', 'CD3G', 'CD8A', 'CD4'],
      'NK': ['EOMES', 'GZMA', 'CCL5', 'XCL1', 'NCAM1'],
      'B': ['CD19', 'MS4A1', 'CD79A', 'CD79B'],
      'PLASMA': ['JCHAIN'],
      'ILC2': ['GATA3', 'IL1RL1', 'IL5', 'IL13'],
      'MYELOID': ['LYZ', 'ITGAM'],
      'NEU': ['S100A8', 'S100A9'],
      'DC': ['ITGAX', 'XCR1', 'CLEC9A', 'CLEC10A'],
      'MON': ['VCAN', 'FN1', 'MAFB'],
      'MAC': ['C1QA', 'MRC1', 'CD68', 'PPARG', 'CEBPB', 'BHLHE41', 'P2RY12', 'TMEM119']}
dp = sc.pl.dotplot(adata, var_names=GD, groupby='leiden_8', cmap='RdPu',
                   standard_scale='var', return_fig=True)
dp.add_totals().show()

# =============================================================================
# 4. save processed object
#    convert to human_sepsis_GSE279452_887650c_processed.rds for the R part
# =============================================================================
adata.write_h5ad(f'human_sepsis_GSE279452_{adata.n_obs}c_processed.h5ad')
