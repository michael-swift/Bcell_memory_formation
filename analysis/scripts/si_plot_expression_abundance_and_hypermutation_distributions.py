import sys, glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
sys.path.insert(0, './analysis/scripts/')
from tb_colormaps import * 

plt.style.use('analysis/scripts/bursa.mplstyle')

######################## PATH CONFIG ################################

df_loc = 'data/vdj/integrated_cell_calls_ambient_annotated.gz'
FIGURE_OUTDIR = 'analysis/figures/si figures/'

######################## ############ ################################

tissue_indices = {'P':0,
                  'B':1,
                  'L':2,
                  'S':3
                 }

fig, ax = plt.subplots(4,3, figsize = (6.5,5))

ax = ax.flatten()
ax_expression = [ax[0], ax[3], ax[6], ax[9]]
ax_abundances = [ax[1], ax[4], ax[7], ax[10]]
ax_mutations = [ax[2], ax[5], ax[8], ax[11]]

df = pd.read_table(df_loc, usecols = ['sample_uid',
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
                                  'is_ambient_source'], 
                                                nrows=None)
df = df[df.sample_uid_gex.notna()]
print('Total cells in GEX samples: ', df.shape[0])
print(df.contaminant_status.value_counts())
df = df[df.contaminant_status.isin(['unlikely_contaminant','source'])]
print('After dropping contaminants and VDJ-free droplets: ', df.shape[0])
df = df[(df['probable_hq_single_b_cell']==True)
        & (df.multiplet_status==1)]
print('Singlet B cells: ', df.shape[0])
acceptable_cells = (~(df.vdj_is_from_ambient==True).astype(bool))
acceptable_cells[df.is_ambient_source] = True
df = df[acceptable_cells]
print('After dropping ambient: ', df.shape[0])


for donor in donor_colors.keys():
    
    print(donor)
    for axis in [ax_mutations]:
        axis[0].plot([],[], lw = 8, label = donor, color = donor_colors[donor])

    called_cell_df = df[df.donor==donor]

    cell_frequencies = called_cell_df.groupby(['tissue', 'vdj_sequence']).size().reset_index()
    cell_frequencies = cell_frequencies.pivot(index='vdj_sequence', columns='tissue', values=0)

    # plot expression and abundance variation
    for _, tissue in enumerate(called_cell_df.tissue.unique()):
        it = tissue_indices[tissue[0]]
        subset = called_cell_df.tissue == tissue

        sns.ecdfplot(data=called_cell_df[subset],
                    x='n_umis',
                    log_scale=(True,True),
                    ax=ax_expression[it],
                    legend=False,
                    complementary=True,
                    color=donor_colors[donor])

        ax_expression[it].set_ylim([10**-5,1.1])

        v_divergence = 100-called_cell_df[subset]['v_pident']
        sns.ecdfplot(data=v_divergence,
                    log_scale=(False,True),
                    legend=False,
                    complementary=True,
                    color=donor_colors[donor],
                    ax=ax_mutations[it])
        ax_mutations[it].axvline(x=v_divergence.mean(), 
                                color=donor_colors[donor], 
                                zorder=10, 
                                lw=1)
        ax_mutations[it].set_ylim([10**-5,1.1])

        sns.ecdfplot(data=cell_frequencies,
                    x=tissue,
                    log_scale=(True,True),
                    legend=False,
                    complementary=True,
                    color=donor_colors[donor],
                    ax=ax_abundances[it])
        ax_abundances[it].set_ylim([10**-5,1.1])

for tissue in ['PB', 'BM', 'SP', 'LN']:
    it = tissue_indices[tissue[0]]
    for axis in [ax_abundances, ax_expression, ax_mutations]:
        axis[it].set_yscale('log', nonpositive='clip')
        axis[it].set_ylabel('')
        axis[it].text(0.9, 0.9, tissue, weight='bold', transform = axis[it].transAxes)
        axis[it].set_xlabel('')
        ax_expression[it].set_ylabel('Proportion greater')

    ax_abundances[it].set_xlim([0.5,200])

ax_abundances[3].set_xlabel('Num. of B cells with VDJ')
ax_expression[3].set_xlabel('IGH UMIs attributable to cell')
ax_mutations[3].set_xlabel('Divergence from germline variant\n(percent)')

fig.legend(loc='upper left', bbox_to_anchor=(1.0,1), frameon=False, handlelength=0.5)
fig.tight_layout()
sns.despine(fig)
fig.savefig(f'{FIGURE_OUTDIR}/basic_vdj_stats_by_tissue.pdf', bbox_inches='tight')

