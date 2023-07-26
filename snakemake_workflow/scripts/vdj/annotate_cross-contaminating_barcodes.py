import sys
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.style.use('./config/bursa.mplstyle')

######################## PATH CONFIG ################################

parser = argparse.ArgumentParser()
parser.add_argument('-input_paths', nargs='+', required="True",
                    help="list of files to combine")
parser.add_argument('-outname', required="True", 
                    help="name for combined sample set")
parser.add_argument('-outdir', default='.')
parser.add_argument('-figure_outdir', default='.')
parser.add_argument('-locus', default='IGH')

args = parser.parse_args()

filenames = args.input_paths
outname = args.outname
outdir = args.outdir
FIGURE_OUTDIR = args.figure_outdir
TARGET_LOCUS = args.locus

seq_col = 'vdj_sequence'
######################## ############ ################################

all_dfs = []
for filename in filenames: 

    df = pd.read_table(filename, nrows=None)
    df = df[df.locus==TARGET_LOCUS]

    all_dfs.append(df)
df = pd.concat(all_dfs)

unique_sample_uids = sorted(df.sample_uid.unique().tolist())

pcc = df.groupby(['sample_uid', 'cb', seq_col]).size().reset_index()
pcc = pcc.pivot(index = ['cb', seq_col], columns='sample_uid', values=0).notna()
num_samples = pcc.sum(axis=1)
single_sample_counts = pcc[num_samples==1]
single_sample_counts = single_sample_counts.sum()
multisample_cell_indices = pcc[num_samples>1].reset_index()
multisample_cell_indices = multisample_cell_indices['cb'] + "+" + multisample_cell_indices[seq_col]
multisample_cell_indices = multisample_cell_indices.to_list()

sys.stderr.write("Number of unique CB-VDJ combinations by the number of unique samples they were found in:\n")
print(num_samples.value_counts(), file=sys.stderr)

####### MAKE A HEATMAP ##########
pcc = pcc.groupby(list(pcc.columns)).size()
pcc = pcc.reset_index().rename(columns={0:'n_cb_vdjs'})
pcc['n_samples'] = pcc[unique_sample_uids].sum(axis=1)
single_sample = pcc[pcc.n_samples == 1]
pcc = pcc[pcc.n_samples==2]

for col in unique_sample_uids:
    if col in pcc:
        pcc.loc[pcc[col]==True, col] = col
        pcc.loc[pcc[col]==False, col] = np.nan
pcc['first'] = pcc.apply(lambda x: x[x.notna()].iloc[0], axis=1)
pcc['second'] = pcc.apply(lambda x: x[x.notna()].iloc[1], axis=1)

pcc = pcc.pivot(index='first', columns='second', values='n_cb_vdjs').fillna(0)
pcc = pcc.reindex(unique_sample_uids, fill_value=0)
for col in unique_sample_uids:
    if col in pcc:
        pass
    else:
        pcc[col] = np.zeros(len(unique_sample_uids))
pcc = pcc[unique_sample_uids]

for col in unique_sample_uids:
    pcc.loc[col,col] = single_sample_counts.loc[col]

fig, ax = plt.subplots(figsize=(8.1,7))
cbar_kws=dict(use_gridspec=False, pad=0.01,shrink=0.5, label=r'$\log_{10}(n_\mathrm{CB-VDJs} + 1)$')

sns.heatmap(np.log(pcc+1)/np.log(10), vmax=5, cmap='YlGn', cbar_kws=cbar_kws)
ax.set_title('Number of likely cross-contaminating cell-associated VDJs')

ax.set_xlabel('')
ax.set_ylabel('')

# Show all ticks and label them with the respective list entries
ax.set_xticks(np.arange(len(pcc.index))+0.5 , labels=pcc.index)
ax.set_yticks(np.arange(len(pcc.columns))+0.5, labels=pcc.columns)

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=90, ha="right", va="center",
         rotation_mode="anchor")

# Loop over data dimensions and create text annotations.
for i in range(len(pcc.index)):
    for j in range(len(pcc.columns)):
        val = int(pcc.iloc[i, j])
        if (val > 0) and (i!=j):
            text = ax.text(j+0.5, i+0.5, val,ha="center", va="center", color="k", fontsize=4)
# PLACE GRIDLINES TO DENOTE DONORS
donor_sample_counts = df.groupby('donor')['sample_uid'].nunique()
donor_sample_counts = donor_sample_counts.cumsum()

for donor in donor_sample_counts:
    ax.axvline(donor, lw=1.5, color='w')
    ax.axhline(donor, lw=1.5, color='w')

fig.savefig(f'{FIGURE_OUTDIR}/{TARGET_LOCUS}_cross_contaminants_heatmap.pdf', bbox_inches='tight')
plt.close(fig)

####### PLOT READ AND UMI STATS ##########

df['cell_id'] = df['cb'] + "+" + df[seq_col]

cdf = df[df.cell_id.isin(multisample_cell_indices)]
ncdf = df[~df.cell_id.isin(multisample_cell_indices)]

cdf = cdf.pivot(index='cell_id', columns='sample_uid', values='n_umis')
cdf['min'] = cdf.min(axis=1)
cdf['max'] = cdf.max(axis=1)
cdf['total'] = cdf['min'] + cdf['max']
cdf['frac'] = cdf['min']/cdf['total']

fig, ax = plt.subplots(1,3, figsize=(6.5,2.2))

counts_r, bins_r = np.histogram(cdf['frac'], bins = np.logspace(-4,0,50))
counts_t, bins_t = np.histogram(cdf['total'], bins = np.logspace(0,5,50))

ax[1].plot(bins_r[1:], (counts_r), lw=1)
ax[0].plot(bins_t[1:], (counts_t), lw=1)
ax[2].scatter(cdf['total'], cdf['frac'], marker='o', alpha=0.2, lw=0)
sns.despine(fig)

for i in [0,1]:
    ax[i].set_xscale('log')
    ax[i].set_ylim([0,40])
for i in [2]:
    ax[i].set_yscale('log')
    ax[i].set_xscale('log')
    ax[i].set_ylim([10**-4,1])

ax[1].set_ylabel('Counts')
ax[1].set_xlabel('Fraction of UMIs in secondary sample')

ax[0].set_ylabel('Counts')
ax[0].set_xlabel('Total UMIs')

ax[2].set_ylabel('Fraction of UMIs in secondary sample')
ax[2].set_xlabel('Total UMIs')

fig.tight_layout()
fig.savefig(f'{FIGURE_OUTDIR}/{TARGET_LOCUS}_cross_contaminants_umi_ratio.pdf', bbox_inches='tight')
plt.close(fig)

contamination_associated_samples = [x for x in cdf.columns if x.startswith('TBd')]
in_contamination_associated_sample_nc = ncdf['sample_uid'].isin(contamination_associated_samples)
ncdf = ncdf[in_contamination_associated_sample_nc]

fig, ax = plt.subplots(1,1,figsize=(4,2))
nc_counts_t, bins_t = np.histogram(ncdf['n_umis'].values, bins = bins_t)
total_c_umis = cdf['min'].sum()
total_nc_umis = cdf['max'].sum() + ncdf['n_umis'].values.sum()

pc = total_c_umis/(total_nc_umis+total_c_umis)
print("Total UMIs associated with likely droplet of origin:", 
      file=sys.stderr)
print(total_nc_umis, file=sys.stderr)
print("Total contaminant:", 
      file=sys.stderr)
print(total_c_umis, file=sys.stderr)
print("Estimated probability of contamination event per UMI:",
      file=sys.stderr)
print(pc, file=sys.stderr)

transfer_probability_per_cell = 1 - np.exp(-bins_t[1:]*pc) - np.exp(-(bins_t[1:]-1)*pc)*bins_t[1:]*pc


ax.plot(bins_t[1:], counts_t, lw=1, label='shared cell barcode-VDJs')
nc_counts_t, bins_t = np.histogram(ncdf['n_umis'].values, bins = bins_t)
ax.plot(bins_t[1:], nc_counts_t, lw=1, label='private cell barcode-VDJs')


ax.set_xlabel('Total IGH UMIs')
ax.set_ylabel('Number of cell barcode-VDJs')

with_theoretical_expectations = True
if with_theoretical_expectations:

    ax.plot(bins_t[1:], transfer_probability_per_cell * nc_counts_t, 
        lw=1, label='model, expected contamination events')


ax.set_yscale('log')
ax.set_xscale('log')
ax.legend(bbox_to_anchor=(0.4,1.05), frameon=False)
ax.set_ylim([0.1, ax.get_ylim()[1]])
fig.tight_layout()
sns.despine(fig)

fig.savefig(f'{FIGURE_OUTDIR}/{TARGET_LOCUS}_cross_contaminants_vs_not_expression_with_pred.pdf', bbox_inches='tight')


########### ANNOTATE CROSS-CONTAMINATING CELLS ######################


cdf = cdf.drop(['min', 'max', 'total'], axis=1)
cdf['source'] = cdf[[x for x in cdf.columns if x in unique_sample_uids]].idxmax(axis=1)
alternative = cdf[[x for x in cdf.columns if x in unique_sample_uids]].idxmin(axis=1)
cdf.loc[(cdf.frac > 10**-2), 
        'source'] = cdf.loc[(cdf.frac > 10**-2), 
                            'source'] + "|" + alternative.loc[cdf.frac>10**-2]
cdf = cdf.reset_index()

cdf = cdf.melt(id_vars=['cell_id','source','frac'], value_name='n_umis')
cdf = cdf[cdf.n_umis.notna()]
cdf = cdf[['cell_id','sample_uid','source']]

cdf['contaminant_status'] = ''
cdf.loc[(cdf.source == cdf.sample_uid),'contaminant_status'] = 'source'
cdf.loc[(cdf.source != cdf.sample_uid),'contaminant_status'] = 'contaminant'
cdf.loc[cdf.source.str.contains("|", regex=False), 'contaminant_status'] = 'unknown_provenance'

cc_cell_dict = cdf.set_index(['cell_id', 'sample_uid']).to_dict()


df['source'] = df.apply(lambda x: 
            cc_cell_dict['source'].get((x.cell_id, x.sample_uid), x.sample_uid), axis=1)
    

df['contaminant_status'] = df.apply(lambda x: 
            cc_cell_dict['contaminant_status'].get((x.cell_id, x.sample_uid), 'unlikely_contaminant'), axis=1)
    
print('Contaminant status assigned:', sys.stderr)
print(df.contaminant_status.value_counts())

df.to_csv(f'{outdir}/{outname}_{TARGET_LOCUS}.tsv.gz', index=False, sep='\t')
