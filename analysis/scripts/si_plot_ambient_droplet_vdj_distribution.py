import sys, glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
sys.path.insert(0, './analysis/scripts/')
from tb_colormaps import * 

plt.style.use('analysis/scripts/bursa.mplstyle')
from scipy.stats import poisson
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
                                    'v_mismatch'],
                                low_memory=False)

################################################################################################################
# Focus on gex samples
df = df[df.sample_uid_gex.notna()]

df['has_vdj'] = df.vdj_sequence.notna()

def gini(array):
    array = np.sort(array) 
    index = np.arange(1,array.shape[0]+1) 
    n = array.shape[0]
    return ((np.sum((2 * index - n  - 1) * array)) / (n * np.sum(array))) #Gini coefficient

total_volume_L = 150 * 10**-6
individual_droplet_volume = 100* 10**-12
n_droplets = total_volume_L/individual_droplet_volume

##### plot droplet distributions within categories

is_ambient_vdj = df[(df.has_vdj)].groupby(['sample_uid_gex', 'vdj_sequence'])['probable_hq_single_not_b_cell'].sum()
is_ambient_vdj = is_ambient_vdj > 0
is_ambient_vdj = is_ambient_vdj.to_dict()

df['vdj_is_ambient'] = False
df.loc[df.has_vdj, 'vdj_is_ambient'] = df[df.has_vdj].apply(lambda x: is_ambient_vdj[(x.sample_uid_gex, x.vdj_sequence)], axis=1)

ambient_vdjs = df[df.vdj_is_ambient].groupby(['sample_uid_gex',
                                              'vdj_sequence'])['n_umis'].agg(list)
ambient_vdjs = ambient_vdjs.map(lambda x: np.asarray(x)).reset_index()
ambient_vdjs['n_droplets'] = ambient_vdjs['n_umis'].map(lambda x: len(x))
ambient_vdjs['total'] = ambient_vdjs['n_umis'].map(lambda x: sum(x))
ambient_vdjs['max'] = ambient_vdjs['n_umis'].map(lambda x : max(x))
ambient_vdjs['max_f'] = ambient_vdjs['n_umis'].map(lambda x : max(x))/ambient_vdjs['total']
ambient_vdjs['mean'] = ambient_vdjs['total']/ambient_vdjs['n_droplets']
ambient_vdjs['median'] = ambient_vdjs['n_umis'].map(lambda x : np.median(x))

ambient_vdjs['gini'] = ambient_vdjs['n_umis'].map(lambda x: gini(x))

fig, axes = plt.subplots(1,3,figsize=(6,2.5), sharey=True)

axes[0].scatter(ambient_vdjs['n_droplets'], ambient_vdjs['total'], 
                marker = '.', 
                lw=0, 
                clip_on=False, 
                color='k', alpha=0.2)

axes[1].scatter(ambient_vdjs['max_f'], ambient_vdjs['total'], 
                marker = '.', 
                lw=0, 
                clip_on=False, 
                color='k', alpha=0.2)
multidroplet = ambient_vdjs['n_droplets'] > 1
axes[2].scatter(ambient_vdjs.loc[multidroplet,'gini'], 
                ambient_vdjs.loc[multidroplet,'total'], 
                marker = '.', 
                lw=0, 
                clip_on=False, 
                color='k', alpha=0.2)
for ax in axes:
    ax.set_yscale('log')
    ax.set_xscale('log')
    xs = np.logspace(0,4,10)
    ax.set_ylim([1,10**6])
axes[0].set_xlim([1,10**6])
axes[1].set_xlim([10**-2,1])
axes[2].set_xlim([10**-2,1])
line1,= axes[0].plot(xs, 2*xs, color='r', lw=0.5, alpha=0.5)
line1.set_dashes((3,2))

ys = xs
line1,= axes[1].plot(2/ys, ys, color='r', lw=0.5, alpha=0.5)
line1.set_dashes((3,2))


axes[0].set_xlabel('Number of droplets associated\nwith ambient VDJ')
axes[1].set_xlabel('Fraction of UMIs in largest droplet\nwith ambient VDJ')
axes[2].set_xlabel('Gini coefficient')

axes[0].set_ylabel('Total UMIs associated\nwith ambient VDJ')
axes[1].set_ylabel('')
axes[2].set_ylabel('')
sns.despine(fig)
fig.tight_layout()
fig.savefig(f'{FIGURE_OUTDIR}/droplet_distribution_of_ambient_vdjs.pdf', bbox_inches='tight')


##### now look into ambient droplets in more detail

df['not_switched'] = df['c_call'].isin(['IGHM', 'IGHD'])

print("unique vdjs associated with non-b:", df[df['vdj_is_ambient']]['vdj_sequence'].nunique())

import itertools

MIN_UMIS_ALL = [0,10,20,50,100,200,500]
FRAC_TOTAL_ALL = [0,0.25,0.4,0.5,0.75,0.8,0.9,0.95,0.99]
search_space = list(itertools.product(MIN_UMIS_ALL, FRAC_TOTAL_ALL))
TARGET_STATS =['n_cells',
             'frac_droplets_accounted',
             'frac_vdjs_accounted',
             'frac_umis_accounted',
             'frac_w_gex',
             'frac_b_lineage',
             'frac_plasma',
             'frac_single_droplet']


possible_source_stats = pd.DataFrame({'min_umis':[x for x,y in search_space],
                                      'frac_total':[y for x,y in search_space]})
for stat in TARGET_STATS:
    possible_source_stats[stat] = np.nan



ambient_df = df[df.vdj_is_ambient].copy()

ambient_vdjs = ambient_vdjs.set_index(['sample_uid_gex', 'vdj_sequence']).to_dict()
ambient_df['has_gex'] = ambient_df.Immune_All_Low_predicted_labels.notna()

ambient_df['vdj_umi_max'] = ambient_df.apply(
                                           lambda x: ambient_vdjs['max'].get(
                                                     (x.sample_uid_gex, x.vdj_sequence)),
                                           axis=1)
ambient_df['vdj_umi_total'] = ambient_df.apply(
                                           lambda x: ambient_vdjs['total'].get(
                                                     (x.sample_uid_gex, x.vdj_sequence)),
                                           axis=1)
ambient_df['vdj_umi_median'] = ambient_df.apply(
                                       lambda x: ambient_vdjs['median'].get(
                                                 (x.sample_uid_gex, x.vdj_sequence)),
                                       axis=1)

ambient_df['larger_than_median'] = (ambient_df.n_umis 
                                    > ambient_df.vdj_umi_median)

ambient_df['is_ASC'] = ambient_df['Immune_All_Low_predicted_labels'].str.startswith("Plasma")
has_ASC_gex = ambient_df.groupby(['sample_uid_gex','vdj_sequence'])['is_ASC'].sum()

print('total ASC = ', has_ASC_gex.sum()) 

has_ASC_gex = has_ASC_gex > 0

ASC_vdjs = has_ASC_gex[has_ASC_gex].reset_index().vdj_sequence.values
print('unique VDJs among ASCs =', len(np.unique(ASC_vdjs)))

print('percent vdjs identifiable ASC = ', has_ASC_gex.sum()/has_ASC_gex.shape[0]) 
print('total sample_uid_gex by vdj combos = ', has_ASC_gex.shape[0])

TOTAL_VDJS = ambient_df.vdj_sequence.nunique()
TOTAL_DROPLETS = ambient_df.shape[0]
print('total droplets:', TOTAL_DROPLETS)
TOTAL_UMIS = ambient_df.n_umis.sum()

has_ASC_vdj = ambient_df['vdj_sequence'].isin(ASC_vdjs)
fraction_of_umis_explainable = ambient_df.loc[has_ASC_vdj, 'n_umis'].sum()/TOTAL_UMIS
fraction_of_droplets_explainable = has_ASC_vdj[has_ASC_vdj].shape[0]/TOTAL_DROPLETS

print('percent droplets from ASC = ', fraction_of_droplets_explainable) 
print('percent UMIs from ASC = ', fraction_of_umis_explainable) 

for it, row in possible_source_stats.iterrows():

    # source_filter = (ambient_df.larger_than_median.astype(bool) 
    #                 & (ambient_df.n_umis > row.min_umis).astype(bool)
    #                 & (ambient_df.n_umis > (row.frac_total * ambient_df.vdj_umi_total)).astype(bool))

    source_filter = ((ambient_df.n_umis == ambient_df.vdj_umi_max).astype(bool) 
                    & (ambient_df.n_umis > row.min_umis).astype(bool)
                    & (ambient_df.n_umis > (row.frac_total * ambient_df.vdj_umi_total)).astype(bool))


    subset = ambient_df[source_filter]
    if subset.shape[0]>0:
        subset_vdjs = subset.vdj_sequence.unique()
        droplets_accounted = sum(1 for x in ambient_df.vdj_sequence.values
                             if x in subset_vdjs)
        vdjs_accounted = sum(1 for x in ambient_df.vdj_sequence.unique()
                             if x in subset_vdjs)
        umis_accounted = ambient_df[ambient_df.vdj_sequence.isin(subset_vdjs)]['n_umis'].sum()
        gex_status_fractions = subset.has_gex.value_counts()/subset.has_gex.value_counts().sum()
        celltype_fractions = (subset['Immune_All_Low_predicted_labels'].value_counts()/
                              subset['Immune_All_Low_predicted_labels'].value_counts().sum())
        droplet_numbers = subset.groupby(['sample_uid_gex', 'vdj_sequence']).size().value_counts()

        frac_w_gex = gex_status_fractions[True]
        bcells = celltype_fractions.index.str.contains("B cell|Plasma|B-")
        plasma = celltype_fractions.index.str.contains("Plasma")
        # print(celltype_fractions[bcells])
        bcell_fraction = celltype_fractions[bcells].sum()
        plasma_fraction = celltype_fractions[plasma].sum()
        single_droplet_fraction = droplet_numbers[1]/droplet_numbers.sum()

        possible_source_stats.at[it,'n_cells'] = subset.shape[0]
        possible_source_stats.at[it, 'frac_droplets_accounted'] = droplets_accounted/TOTAL_DROPLETS
        possible_source_stats.at[it, 'frac_vdjs_accounted'] = vdjs_accounted/TOTAL_VDJS        
        possible_source_stats.at[it, 'frac_umis_accounted'] = umis_accounted/TOTAL_UMIS
        possible_source_stats.at[it, 'frac_w_gex'] = frac_w_gex
        possible_source_stats.at[it, 'frac_b_lineage'] = bcell_fraction
        possible_source_stats.at[it, 'frac_plasma'] = plasma_fraction
        possible_source_stats.at[it, 'frac_single_droplet'] = single_droplet_fraction
    else:
        possible_source_stats.at[it,'n_cells'] = 0
        possible_source_stats.at[it, 'frac_droplets_accounted'] = np.nan 
        possible_source_stats.at[it, 'frac_vdjs_accounted'] = np.nan
        possible_source_stats.at[it, 'frac_umis_accounted'] = np.nan
        possible_source_stats.at[it, 'frac_w_gex'] = np.nan
        possible_source_stats.at[it, 'frac_b_lineage'] = np.nan
        possible_source_stats.at[it, 'frac_plasma'] = np.nan
        possible_source_stats.at[it, 'frac_single_droplet'] = np.nan

possible_source_stats.to_csv('data/vdj/ambient_droplet_stats_max_only.tsv', sep='\t')

possible_source_stats = pd.read_table('data/vdj/ambient_droplet_stats_max_only.tsv')
fig, ax = plt.subplots(2,4, figsize=(int(2*len(FRAC_TOTAL_ALL)),
                                     int(2*len(MIN_UMIS_ALL)/2)))
ax = ax.flatten()

for it, x in enumerate(TARGET_STATS):
    if x.startswith("n"):
        fmt='.4g'
    else:
        fmt='.2g'
    data = possible_source_stats.pivot(columns='frac_total',
                                       index='min_umis',
                                       values=x)
    sns.heatmap(data, cmap='YlGn', ax = ax[it], annot=True, fmt=fmt)
    ax[it].set_title(x)
fig.tight_layout()
fig.savefig(f'{FIGURE_OUTDIR}/potential_source_ambient_droplet_stats_max_only.pdf', bbox_inches='tight')


## First compare statistics of ambient VDJs with and without identifiable plasma cell sources

print('Comparing statistics of ambient VDJs with identifiable plasma cell sources')
print('and statistics of ASC VDJs not found in non-B cell')

asc_associated = ambient_df.apply(lambda x: has_ASC_gex[(x.sample_uid_gex, x.vdj_sequence)] , axis=1)
n_droplets_with_ASC = ambient_df[asc_associated& df['probable_hq_single_b_cell']].groupby(['sample_uid_gex','vdj_sequence'])['n_umis'].agg(list).map(lambda x: np.asarray(x))

is_ASC =df.Immune_All_Low_predicted_labels.str.contains("Plasma")
ASC_vdj_list = df[is_ASC==True].groupby(['sample_uid_gex', 'vdj_sequence'])['n_umis'].agg(list).map(lambda x: np.asarray(x))
ASC_vdj_list = ASC_vdj_list.reset_index().groupby(['sample_uid_gex']).agg(list)
ASC_vdj_list = ASC_vdj_list['vdj_sequence'].to_dict()
has_ASC_vdj = df.apply(lambda x: x.vdj_sequence in ASC_vdj_list.get(x.sample_uid_gex, [])
                       ,axis=1)
ASC_filter = has_ASC_vdj & ~(df.vdj_is_ambient.astype(bool)) & df['probable_hq_single_b_cell']

# Collect per-vdj stats

all_ASC_vdjs = df[ASC_filter]
n_droplets_other_ASC_VDJS = all_ASC_vdjs.groupby(['sample_uid_gex','vdj_sequence'])['n_umis'].agg(list).map(lambda x: np.asarray(x))



fig, axes = plt.subplots(1,2,figsize=(6.5,1.8))

for item, color, label in [(n_droplets_other_ASC_VDJS,'k', 'ASC but not non-B cell'),
                           (n_droplets_with_ASC, 'red', 'ASC and non-B cell'), 
                    ]:
    item = item.reset_index()
    item['n_droplets'] = item['n_umis'].str.len()
    item['total'] = item['n_umis'].map(lambda x: np.sum(x))

    item['max_f'] = item['n_umis'].map(lambda x: np.max(x))/item['n_umis'].map(lambda x: np.sum(x))
    item['gini'] = item['n_umis'].map(lambda x: gini(x))
    sns.ecdfplot(item, x='n_droplets', ax =axes[0], color=color, complementary=True)
    sns.ecdfplot(item[(item.n_droplets>1)], x='gini', complementary=False, ax=axes[1], color=color)
    axes[1].plot([],[],lw=5, label=label, color=color)
axes[1].legend(frameon=False, bbox_to_anchor=(1,1))

axes[0].set_ylabel('Proportion greater')
axes[0].set_xlabel('Number of droplets with VDJ')
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[0].set_ylim([axes[0].get_ylim()[0], 1.1])
axes[1].set_xscale('log')
axes[1].set_yscale('log')
axes[1].set_ylabel('Proportion smaller')
axes[1].set_xlabel('Gini coefficient')

fig.tight_layout()
sns.despine(fig)
fig.savefig(f'{FIGURE_OUTDIR}/droplet_distribution_per_ASC_VDJ_class.pdf', bbox_inches='tight')

