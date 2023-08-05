import sys, glob

import numpy as np
import scipy 
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
from tb_colormaps import * 
plt.style.use('analysis/scripts/bursa.mplstyle')

sys.path.insert(0, './analysis/scripts/')

######################## PATH CONFIG ################################

df_loc = 'data/vdj/integrated_cell_calls_ambient_annotated.tsv.gz'
ASC_subtype_df_loc = 'data/vdj/ASC_subtypes.tab'
FIGURE_OUTDIR = 'analysis/figures/si figures/'

######################## ############ ################################


df = pd.read_table(df_loc, usecols = ['cb',
                                  'sample_uid',
                                  'sample_uid_gex',
                                  'vdj_sequence',
                                  'n_umis',
                                  'v_pident',
                                  'donor',
                                  'tissue',
                                  'contaminant_status',
                                  'multiplet_status',
                                  'probable_hq_single_b_cell',
                                  'vdj_is_from_ambient',
                                  'is_ambient_source',
                                  'Immune_All_Low_predicted_labels'], 
                                                nrows=None)
df = df[df.sample_uid_gex.notna()]
print('Total cells in GEX samples: ', df.shape[0])
df = df[(df['probable_hq_single_b_cell']==True)
        & (df.multiplet_status==1)]
print('Singlet B cells: ', df.shape[0])
acceptable_cells = (~(df.vdj_is_from_ambient==True).astype(bool))
acceptable_cells[df.is_ambient_source] = True
df = df[acceptable_cells]
print('After dropping ambient: ', df.shape[0])
df = df[df.donor=='TBd3']
print('Focusing on TBd3: ', df.shape[0])

df['sample_preparation'] = df.sample_uid.map(lambda x: x.split("_")[1])

print(df.sample_preparation.value_counts())
#add subtype annotations

celltypist_simpler = {
    "Proliferative germinal center B cells": "GC B cells",
    "Germinal center B cells": "GC B cells", 
    "Age-associated B cells" : "ABCs"
}

ASC_subtype_df = pd.read_table(ASC_subtype_df_loc)

ASC_subtype_df.columns=['barcode','subtype','sample_uid']
ASC_subtype_df['cb'] = ASC_subtype_df.barcode.str.split("-").map(lambda x: x[0])
ASC_subtype_df['cb_suid'] = ASC_subtype_df['cb'] + "_" + ASC_subtype_df['sample_uid']
ASC_subtype_dict = ASC_subtype_df.set_index('cb_suid')['subtype'].to_dict()
df['cb_suid'] = df['cb'] + "_" + df['sample_uid']
df['ASC_subtype'] = df.cb_suid.map(ASC_subtype_dict)
df['celltype'] = df['Immune_All_Low_predicted_labels'].copy()
df['celltype'] = df.celltype.map(lambda x: celltypist_simpler.get(x,x))
df.loc[df.ASC_subtype.notna(), 'celltype'] = df['ASC_subtype']
common_celltypes = df['celltype'].value_counts()
common_celltypes = common_celltypes.index
# drop low-quality ASCs that did not get assigned a subtype
common_celltypes = [x for x in common_celltypes if not x.startswith('Plasma')]

total_counts = df[df.celltype.isin(common_celltypes)].groupby(['donor','sample_preparation','tissue']).size()
celltype_dist = df[df.celltype.isin(common_celltypes)].groupby(['donor','sample_preparation','tissue'])[['celltype']].value_counts(normalize=True).reset_index()
celltype_dist = celltype_dist.set_index(['donor', 'sample_preparation', 'tissue']).rename(columns={0:'fraction'})
celltype_dist['n'] = celltype_dist.index.map(total_counts)
celltype_dist = celltype_dist.reset_index()

CI = scipy.stats.binom.interval(0.95, n= celltype_dist.n, p = celltype_dist.fraction)

celltype_dist['95_CI_lower'] = CI[0]/celltype_dist['n']
celltype_dist['95_CI_upper'] = CI[1]/celltype_dist['n']
celltype_dist['fraction'] = celltype_dist['fraction'] + np.maximum(1/celltype_dist['n'], 10**-4)


fig, ax = plt.subplots(1,1, figsize=(2.7,1.8), sharex=True)

donor3 = celltype_dist[celltype_dist.donor=='TBd3']

yvals = {'Naive B cells':0, 
         'Memory B cells':1,
         'ABCs':2,
         'ASC-1':3,
         'ASC-2':4,
         'ASC-3':5,
         'ASC-4':6,
         'Plasma cells':7,
         'Plasmablasts':8}
yticks = sorted(yvals.values())
yval_lookup = {v:k for k, v in yvals.items()}
yticklabels = [yval_lookup[y] for y in yticks]
ax.set_yticks(ticks=yticks, labels=yticklabels)

dy = {'BM':0.12, 
      'PB':-0.12,
      'fresh':-0.06,
      'frozen':0.06}
kwargs = {'fresh':dict(marker='o', s=10),
          'frozen':dict(marker=(4,1,0), s=18)}

for it, row in donor3.iterrows():
    y = yvals[row.celltype] + dy[row.tissue] + dy[row.sample_preparation]
    color = tissue_colors[row.tissue]
    symbol_kwargs = kwargs[row.sample_preparation]

    ax.plot([row['95_CI_lower'], row['95_CI_upper']],
            [y,y], 
            color=color,
            lw=1, 
            alpha=0.8)
    ax.scatter(row['fraction'],
               y, 
               color=color, 
               **symbol_kwargs, 
               linewidths=0)
    
ax.plot([],[],lw=5, color=tissue_colors['PB'], label='PB')
ax.plot([],[],lw=5, color=tissue_colors['BM'], label='BM')
ax.scatter([],[],  lw=0, color='k', label='fresh', **kwargs['fresh'],)
ax.scatter([],[],  lw=0, color='k', label='frozen', **kwargs['frozen'],)

ax.legend(bbox_to_anchor=(1,1), frameon=False, handlelength=0.5)
ax.set_xscale('log')
ax.set_xlim([10**-4,1])
ax.set_xlabel('Fraction of cells')
ax.set_ylabel('B cell type')
fig.tight_layout()
sns.despine(fig)

fig.savefig(f'{FIGURE_OUTDIR}/freezing_bias.pdf', bbox_inches='tight')
