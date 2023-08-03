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
FIGURE_OUTDIR = 'analysis/si figures/'

######################## ############ ################################

df = pd.read_table(df_loc, usecols=['sample_uid',
                                    'tissue',
                                    'n_umis',
                                    'n_reads',
                                    'vdj_sequence'])
df=df[df.vdj_sequence.notna()]
tissue_dict = df[['sample_uid','tissue']].drop_duplicates()
tissue_dict = tissue_dict.set_index(['sample_uid'])['tissue'].to_dict()

fig, ax = plt.subplots(2, 1, figsize=(3, 4), sharex=True)
# ax_read_counts = ax[0]
ax_umi_counts = ax[1]
ax_umi_quantiles = ax[0]
variable = 'n_umis'

for groupby_var, alpha, lw in [('tissue',1,1), ('sample_uid',0.2, 0.5)]:
    counts = df.groupby([groupby_var,variable]).size()
    counts = counts.reset_index().rename(columns={0:'num'})
    counts = counts.pivot(index=variable, columns=groupby_var, values='num')

    for sample_it, sample in enumerate(counts.columns):
        if groupby_var == 'tissue':
            color = tissue_colors[sample]#tissue_dict[sample]]
        else:
            color = tissue_colors[tissue_dict[sample]]
        kwargs = dict(color = color, lw=lw, alpha=alpha)
        cat_total  = counts[sample].index*counts[sample]
        cat_fraction = cat_total/cat_total.sum()
        inverse_cumsum = np.cumsum(cat_fraction[::-1])[::-1]
        quantile = (counts[sample]/counts[sample].sum(axis=0)).cumsum()

        ax_umi_quantiles.plot(quantile, 1- inverse_cumsum, **kwargs)
        ax_umi_quantiles.set_xlabel(f'Cell quantile (UMIs)')
        ax_umi_quantiles.set_ylabel('Fraction of sample\navailable to smaller cells')
        ax_umi_quantiles.set_yscale('log')
        ax_umi_quantiles.set_ylim([10**-6, 1])
        # ax_umi_quantiles.set_xscale('log')
        ax_umi_quantiles.axvline(x=0.5, color='0.5', lw=1, alpha=0.5)

        count_values = counts.index[counts[sample].notna()]
        quantile_values = quantile[counts[sample].notna()]
        ax_umi_counts.plot(quantile_values, count_values, **kwargs)
        ax_umi_counts.set_xlabel(f'Cell quantile (UMIs)')
        ax_umi_counts.set_ylabel('Number of UMIs')
        ax_umi_counts.set_yscale('log')
        # ax_umi_counts.set_ylim([10**-6, 1])
        # ax_umi_counts.set_xscale('log')
        ax_umi_counts.axvline(x=0.5, color='k', lw=0.5, alpha=0.5)
for tissue in ['PB', 'BM', 'SP', 'LN']:
    ax_umi_quantiles.plot([],[], lw=8, color=tissue_colors[tissue], label=tissue)

ax_umi_quantiles.legend(bbox_to_anchor=(1,1), handlelength=0.6, frameon=False, title='Tissue')
fig.tight_layout()
sns.despine(fig)
fig.savefig('analysis/figures/si figures/read_and_umi_availability.pdf', bbox_inches='tight')
# plt.show()
