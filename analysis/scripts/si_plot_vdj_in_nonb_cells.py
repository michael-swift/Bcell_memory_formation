import sys, glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sys.path.insert(0, './analysis/scripts/')
from tb_colormaps import * 

plt.style.use('analysis/scripts/bursa.mplstyle')

######################## PATH CONFIG ################################

df_loc = 'data/vdj/integrated_cell_calls.tsv.gz'
FIGURE_OUTDIR = 'analysis/figures/si figures/'

######################## ############ ################################

df = pd.read_table(df_loc, usecols=['sample_uid',
                                    'sample_uid_gex',
                                    'sample_type',
                                    'tissue',
                                    'n_umis',
                                    'n_reads',
                                    'vdj_sequence',
                                    'Immune_All_Low_predicted_labels',
                                    'Immune_All_High_predicted_labels',
                                    'probable_hq_single_b_cell',
                                    'possible_b_cell',
                                    'probable_hq_single_not_b_cell',
                                    'rare_or_bad_q_cell',
                                    'multiplet_status'])


################################################################################################################

print(df.sample_type.value_counts())
df['sample_uid'] = df['sample_uid_gex']
df['has_vdj'] = df.vdj_sequence.notna()

tissue_dict = df[['sample_uid','tissue']].drop_duplicates()
tissue_dict = tissue_dict.set_index(['sample_uid'])['tissue'].to_dict()

non_b_cells = df[(df.probable_hq_single_not_b_cell == True)]
plasma_subsets = df.Immune_All_High_predicted_labels.fillna("").str.contains('Plasma')

mean_multiplicity = df[df.multiplet_status>0].groupby('sample_uid')['multiplet_status'].mean()
mean_multiplicity = mean_multiplicity.to_dict()
singlet_rate = df[df.multiplet_status==1].groupby('sample_uid').size()
singlet_rate = singlet_rate/df[df.multiplet_status>0].groupby('sample_uid').size()
singlet_rate = singlet_rate.to_dict()



total_droplets = df[df.Immune_All_High_predicted_labels.notna()]['sample_uid'].value_counts()
fraction_non_b = non_b_cells['sample_uid'].value_counts()/total_droplets
# total_b = (1-fraction_non_b)*total_droplets
total_vdjs = df[df.has_vdj]['sample_uid'].value_counts()
total_vdjs = total_vdjs.to_dict()

total_large_umi_vdjs = df.groupby(['sample_uid','vdj_sequence'])['n_umis'].sum().reset_index()
total_large_umi_vdjs = total_large_umi_vdjs[total_large_umi_vdjs['n_umis']>20].groupby(['sample_uid'])['vdj_sequence'].size()
total_large_umi_vdjs = total_large_umi_vdjs.to_dict()

vdj_rate_in_ambient = non_b_cells.groupby(['has_vdj', 'sample_uid']).size()
vdj_rate_in_ambient = vdj_rate_in_ambient.reset_index()

vdj_rate_in_ambient = vdj_rate_in_ambient.pivot(index='sample_uid', 
                                  columns='has_vdj', 
                                  values=0).fillna(0)

vdj_rate_in_ambient = (vdj_rate_in_ambient.T/vdj_rate_in_ambient.sum(axis=1)).T
vdj_rate_in_ambient['total_vdjs'] = vdj_rate_in_ambient.index.map(total_vdjs)
vdj_rate_in_ambient['tissue'] = vdj_rate_in_ambient.index.map(tissue_dict)
vdj_rate_in_ambient['total_large_umi_vdjs'] = vdj_rate_in_ambient.index.map(total_large_umi_vdjs)

fig, ax = plt.subplots(1,2,figsize=(7,3), sharey=True)
ax = ax.flatten()

vdj_rate_in_ambient['fraction_with_ambient_vdj'] = vdj_rate_in_ambient[True] + 0.5*10 **-3

for it, x, xlabel in [(0, 'sample_uid', 'Unique sample identifier'),
              (1, 'total_vdjs', 'Total number of VDJs with >20 UMIs')]:
    for sample_type, symbol, s in [('fresh','o',20), ('frozen',(4,1,0),30)]:
        subset = vdj_rate_in_ambient.index.str.contains(sample_type)
        g = sns.scatterplot(data=vdj_rate_in_ambient[subset], 
                        y='fraction_with_ambient_vdj', 
                        x=x,
                        hue='tissue',
                        palette=tissue_colors,
                        legend=False,
                        marker=symbol,
                        s=s,
                        ax=ax[it],
                        alpha=0.8,
                        clip_on=False)
        g.set_xlabel(xlabel)
    if x == 'sample_uid':
        g.set_xticklabels(g.get_xticklabels(), rotation=90)
    else:
        g.set_xscale('log')
    g.set_ylabel('Fraction of non-B cell droplets\nwith ambient IGH amplicon')

    g.set_yscale('log')
    g.set_ylim([.5*10**-3, 0.2])

for tissue in ['PB', 'BM', 'LN', 'SP']:
    g.plot([],[], lw=8, color = tissue_colors[tissue], label=tissue)
for sample_type, symbol, s in [('fresh','o',20), ('frozen',(4,1,0),30)]:
    g.scatter([],[], color='k', marker=symbol, s = s, label=sample_type, lw=0)
g.legend(bbox_to_anchor=(1,1), frameon=False, handlelength=0.6)
fig.tight_layout()
sns.despine(fig)
fig.savefig(f'{FIGURE_OUTDIR}/ambient_rates.pdf', bbox_inches='tight')
plt.close(fig)

