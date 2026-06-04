#################################################################
## python code for generating figures related to Figure1.
## Fig1b, s1d, s1e
#################################################################

## load data
import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scanpy as sc
import scipy.sparse as sp
from collections import Counter
import cosg as cosg
import importlib
importlib.reload(cosg)

sc.settings.verbosity = 3
adata = sc.read_h5ad('all/data/mouse_v9_final_120649c.h5ad') 
annot_cat = ['B_pro','B_pre','B_mature','plasma','T_pre','T_naive','T_mature','T_gd','ILC','NK','Neu','Basophil','MastCell','Eos','Mon','Mac','KupfferCell','Langerhans','cDC2','cDC1','pDC','mDC' ]
adata.obs.annot_d = adata.obs.annot_d.cat.reorder_categories(annot_cat)

#################################################################
## Fig. 1b, UMAP plot
del adata.uns['annot_d_colors'] 
col_22 = ['#00417d',
'#ffad30',
'#006d42',
'#7399a0',
'#29a4d4',
'#4e347e',
'#f3762b',
'#956e9e',
'#e1352a',
'#7d9e61',
'#0089ac',
'#f0d969',
'#003e59',
'#fbb396',
'#b22524',
'#8b5f65',
'#e67f98',
'#0065a9',
'#d73659',
'#fa7d56',
'#5a70a8',
'#5ba77a'
]
adata.uns['annot_d_colors'] = col_22
sc.pl.umap(adata, color = 'annot_d', legend_loc=None, title='')

#################################################################
## Fig. s1d, UMAP plot
sc.pl.umap(adata, color= 'annot_r', frameon=False, title = '', legend_loc=None)

#################################################################
## Fig. s1e, heatmap
cosg.cosg(adata,
    key_added='cosg',
        mu=1,
        n_genes_user=50,
               groupby='annot_d')

df_tmp=pd.DataFrame(adata.uns['cosg']['names'][:5,]).T
df_tmp.reindex(annot_cat)
marker_genes_list=np.ravel(df_tmp.reindex(annot_cat))
sc.pl.heatmap(adata, marker_genes_list, swap_axes = True,
             groupby='annot_d',                           
             standard_scale='var',
             cmap='RdPu',save=True,figsize=(20,10), log=True)

