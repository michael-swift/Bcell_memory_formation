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

df_loc = 'data/vdj/integrated_cell_calls.tsv.gz'
FIGURE_OUTDIR = 'analysis/figures/si figures/'

######################## ############ ################################

df = pd.read_table(df_loc, usecols=['sample_uid_gex',
                                    'sample_type',
                                    'tissue',
                                    'n_umis',
                                    'n_reads',
                                    'vdj_sequence',
                                    'Immune_All_Low_predicted_labels',
                                    'probable_hq_single_b_cell',
                                    'probable_hq_single_not_b_cell',
                                    'multiplet_status',
                                    'c_call',
                                    'v_mismatch'])


################################################################################################################

df['has_vdj'] = df.vdj_sequence.notna()

s = {"IGHM":"IGHM|D", "IGHD":"IGHM|D"}
df['switched'] = df.c_call.map(lambda x: s.get(x, 'switched'))
df['switched'] = df['switched'] == 'switched'
df['hypermutated'] = df.v_mismatch.map(lambda x: True if x > 2 else False)
# focus on GEX samples
df = df[df.sample_uid_gex.notna()]

ambient_rate = df[(df.probable_hq_single_not_b_cell == True)].groupby(['sample_uid_gex'])['has_vdj'].mean()
ambient_rate = ambient_rate.to_dict()

########### HQ B cells only ###############

bcells = df[(df.probable_hq_single_b_cell==True) & (df.multiplet_status<2)]
subtype_dict = {'Plasma cells':'ASCs',
                'Plasmablasts':'ASCs',
                'Age-associated B cells':'ABCs',
                'Germinal center B cells':"GC B cells",
                'Proliferative germinal center B cells':'GC B cells'}
df['subtype'] = ''
bcells.loc[:,'subtype'] = bcells.Immune_All_Low_predicted_labels.map(lambda x: subtype_dict.get(x,x))

common_celltypes = bcells['subtype'].value_counts()
common_celltypes = common_celltypes[common_celltypes > 99].index.values

bcells_common = bcells[bcells.subtype.isin(common_celltypes)]
bcells_rare = bcells[~bcells.subtype.isin(common_celltypes).astype(bool)]

for subset, label in [(bcells_common, 'common-singlets')]:
    vdj_rates = subset.groupby(['sample_uid_gex', 'subtype'])['has_vdj'].mean().reset_index()

    # estimate ambient rates only for subtypes represented by at least 10 cells in a sample
    subset_sizes = subset.groupby(['sample_uid_gex','subtype']).size()
    vdj_rates['subset_size'] = vdj_rates.apply(lambda x : subset_sizes[(x.sample_uid_gex, x.subtype)], axis=1)
    vdj_rates = vdj_rates[vdj_rates['subset_size'] > 9]
    
    #further exclude samples if we do not have an independent estimate of the ambient rate
    have_ambient_estimate = vdj_rates.sample_uid_gex.isin(ambient_rate.keys())
    vdj_rates = vdj_rates[have_ambient_estimate]


    ########### B cells odds ratio ###############

    vdj_rates['ambient_rate'] = vdj_rates.sample_uid_gex.map(ambient_rate)
    vdj_rates['p_real_vdj'] = 1 - vdj_rates['ambient_rate']/vdj_rates['has_vdj']
    vdj_rates['p_ambient_vdj'] = vdj_rates['ambient_rate']/vdj_rates['has_vdj'] + 0.5*10**-3


    fig = plt.figure(figsize=(6.5,1.8))
    gs = gridspec.GridSpec(1, 3, figure=fig, width_ratios = [1,1.2,1], wspace=0.4)

    ax_p_any = fig.add_subplot(gs[0])
    ax_p_correct = fig.add_subplot(gs[1], sharey=ax_p_any)
    ax_naive = fig.add_subplot(gs[2])
    g= sns.ecdfplot(data=vdj_rates, 
                    x ='has_vdj', 
                    ax=ax_p_any,
                    hue='subtype',
                    palette=bcelltype_colors_alt,
                    legend=False)

    ax_p_any.set_xlabel("Percent cells in category with VDJ")
    ax_p_any.set_ylabel("Proportion of samples smaller")


    g= sns.ecdfplot(data=vdj_rates, 
                    x ='p_ambient_vdj', 
                    ax = ax_p_correct,
                    hue='subtype',
                    palette=bcelltype_colors_alt, 
                    complementary=True,
                    legend=True)
    ax_p_correct.set_xlim([0.5*10**-3,1.0])
    ax_p_correct.set_xscale('log')
    ax_p_correct.set_xlabel("Probability VDJ is ambient")
    ax_p_correct.set_ylabel("Proportion of samples greater")
    sns.move_legend(g, frameon=False, loc='upper right', bbox_to_anchor=(1.,1.1), title='')

    fig.tight_layout()
    sns.despine(fig)



    ## Naive ambient estimation

    print(vdj_rates.columns)

    for count, marker in [('switched','.')]:#, ('hypermutated', '.')]:# 'hypermutated']:

        data = bcells[(bcells.subtype == "Naive B cells") & (bcells.has_vdj)].groupby("sample_uid_gex")[count].mean()
        data.name = 'measured'
        ambient_baseline = df[(df.probable_hq_single_not_b_cell==True) & (df.has_vdj)].groupby("sample_uid_gex")[count].mean()

        ambient_baseline.name=f'baseline {count}'
        data = pd.merge(data, vdj_rates, 
                        left_on='sample_uid_gex',
                        right_on='sample_uid_gex')
        data = pd.merge(data, ambient_baseline, on='sample_uid_gex')

        data = data[data.subtype == 'Naive B cells']
        data['p_ambient_vdj'] = 1. - data['p_real_vdj']
        data['Ambient-based estimate'] = data.p_ambient_vdj + 10**-4
        data['IGHC-based estimate'] = (data['measured'])/(data[f'baseline {count}']) + 10**-4


        # ax_naive.set_title('Naive B cell contamination rate')
        sns.scatterplot(data,
                        x='Ambient-based estimate', 
                        y='IGHC-based estimate', 
                        color=bcelltype_colors_alt['Naive B cells'], 
                        legend=False, 
                        ax=ax_naive, 
                        marker=marker,
                        clip_on=False)
        ax_naive.set_xlabel('Ambient-based estimate\nof Naive B cell contamination')
        ax_naive.set_ylabel('IGHC-based estimate\nof Naive B cell contamination')

        ax_naive.set_xscale('log')
        ax_naive.set_yscale('log')
        ax_naive.set
        xs = np.logspace(-3.5,-0.5,20,base=10)
        ys_lower=0.2*xs
        ys_upper = 5*xs
        ax_naive.plot(xs,ys_lower,  '--', color='0.5', lw=0.5)
        ax_naive.plot(xs,ys_upper, '--', color='0.5', lw=0.5)
        ax_naive.plot(xs,xs, color='0.5', lw=0.75)
        ax_naive.set_xlim([10**-4,10**0.5])
        ax_naive.set_ylim([10**-4,10**0.5])
        fig.savefig(f'{FIGURE_OUTDIR}/bcell_VDJs_{label}.pdf', bbox_inches='tight')



    vdj_rates.to_csv(f'data/vdj/ambient_rate_and_prob_own_vdj_by_celltype_{label}.tsv',
                     index=True,
                     sep='\t')
