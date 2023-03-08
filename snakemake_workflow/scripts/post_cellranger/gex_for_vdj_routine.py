#!/usr/bin/env python
# coding: utf-8
import glob
import gc
import pandas as pd
import numpy as np
import scanpy as sc
from matplotlib import pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pathlib
import celltypist
import scrublet as scr
from celltypist import models
params = {
    'font.size': 12,
    'axes.titlesize': 12,
    'axes.labelsize': 12,
    'legend.fontsize': 12,
    'xtick.labelsize': 8,
    'ytick.labelsize': 10,
    'font.family': "Arial",
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'figure.dpi': 100
   }
mpl.rcParams.update(params)
sns.set_style("ticks")
sns.set_context(context='paper')
savefig_args = {"dpi": 300, "bbox_inches": "tight", "pad_inches": 0, "transparent": True}
mpl.rc('savefig', dpi=300)



pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 20) 
pd.set_option('display.width', 100)

# Set UP I/O

fig_output_dir = snakemake.params.fig_output
gex_labels_output = snakemake.output.gex_labels
h5ad_outputs = snakemake.output.h5ads
h5ad = snakemake.input
tissue_labels = snakemake.parmas.tissue_labels

def setup_figure_outs(tissue):
    output_dir='{}/figures/QCandAnnotation/gex_{}'.format(output_dir, tissue)
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_suffix = ""
    output_formats = [".png", ".svg"]
    sc.settings.figdir = output_dir
    sc.set_figure_params(format='pdf', transparent=True,)

    return output_dir, output_formats, output_suffix

def run_gex_label_routine(data_path, tissues):
    for tissue in tissues:
        output_dir, output_formats, output_suffix = setup_figure_outs(tissue=tissue)
        def save_figure(fig, name, output_dir=output_dir, output_suffix=output_suffix, output_formats=output_formats, savefig_args=savefig_args):
            for output_format in output_formats:
                fig.savefig(output_dir + "/" + name + output_suffix + output_format, **savefig_args)
        if tissue in file:
                adata = sc.read_h5ad(file)
                # remove IGH and IGL variable genes from highly variable genes for clustering analysis 
                adata.var.loc[adata.var.index.str.contains("IGH|IGL|IGK|FOS|JUN|HSP|RPL"), 'highly_variable'] = False
                # ad hoc exclusion of weird samples
                adata = adata[adata.obs.sample_uid != 'TBd3_fresh_B200']
                adata = cluster(adata, batch_correct=False)
                filter_low_abundance_cell_groups = False
                cell_group = "predicted_labels"
                if filter_low_abundance_cell_groups:
                    select = adata.obs[cell_group].value_counts() > (adata.obs.shape[0] / 1000)
                    adata = adata[adata.obs[cell_group].isin(select[select == True].index)]
                # create low resolution leiden for majority voting
                # High Resolution
                predictions = celltypist.annotate(adata, model = 'Immune_All_Low.pkl', majority_voting=False)
                adata = predictions.to_adata(prefix="Immune_All_Low_")
                # not sure why this isn't added automatically
                adata.uns['log1p'] = {"base":np.e}
                # Low Resolution
                predictions = celltypist.annotate(adata, model = 'Immune_All_High.pkl', majority_voting=False)
                adata = predictions.to_adata(prefix="Immune_All_High_")
                # Majority Voting with my Leiden
                sc.tl.leiden(adata, resolution=0.2)
                predictions = celltypist.annotate(adata, model = 'Immune_All_Low.pkl', majority_voting=True, over_clustering=adata.obs.leiden)
                adata = predictions.to_adata(prefix="my_leiden_")
                # not sure why this isn't added automatically
                adata.uns['log1p'] = {"base":np.e}
                g = sns.jointplot(data = adata.obs, x = "Immune_All_Low_conf_score", y = "Immune_All_High_conf_score", kind = 'kde')
                save_figure(g.figure, "confidence_scores_celltypist")
                # doublet removal and scoring
                for layer in ['umi_counts', 'background_removed']: 
                    scrub = scr.Scrublet(adata.to_df(layer=layer))
                    prediction = 'predicted_doublets_{}'.format(layer)
                    score = 'doublet_scores_{}'.format(layer)
                    doublet_scores, predicted_doublets = scrub.scrub_doublets()
                    adata.obs.loc[:, score] = doublet_scores 
                    adata.obs.loc[:,prediction] = predicted_doublets
                    # modify predicted doublets based on manual score cutoff
                    sns.displot(data = adata.obs, x = score)
                    plt.yscale('log')
                    sns.displot(data = adata.obs, x = score, kind = 'ecdf')
                    adata.obs.loc[:,prediction] = adata.obs[score] > 0.2
                    adata.obs.loc[:,prediction] = adata.obs[prediction].astype(str)
                ## Technical UMAPs
                variables = ['predicted_doublets_umi_counts', 'predicted_doublets_background_removed',
                             'doublet_scores_umi_counts', 'doublet_scores_background_removed', 'Immune_All_High_predicted_labels', 
                             'Immune_All_Low_predicted_labels', 'leiden', 'sample_uid']
                for var in variables:
                    sc.pl.umap(adata, color = var, size = 10, save = "{}_{}".format(var, tissue))

                g = sns.heatmap(sc.metrics.confusion_matrix('predicted_doublets_umi_counts', 'predicted_doublets_background_removed', adata.obs, normalize=False), annot=True)
                save_figure(g.figure, "confusion_matrix")
                    
                ## Plot Background Removal Examples
                gene = "IGKC"
                _umi = sc.get.obs_df(adata, keys = gene, layer='umi_counts')
                _bg = sc.get.obs_df(adata, keys = gene, layer='background_removed')
                data = pd.concat([_bg, _umi], axis=1)
                data.columns = ['cbender', 'umi']
                data = data.melt()
                g = sns.displot(data, x = 'value', hue = 'variable', kind = 'ecdf', log_scale=True)
                save_figure(g.figure, "{}_bg_remove".format(gene))

                gene = "AZU1"
                _umi = sc.get.obs_df(adata, keys = gene, layer='umi_counts')
                _bg = sc.get.obs_df(adata, keys = gene, layer='background_removed')
                data = pd.concat([_bg, _umi], axis=1)
                data.columns = ['cbender', 'umi']
                data = data.melt()
                g = sns.displot(data, x = 'value', hue = 'variable', kind = 'ecdf', log_scale=True)
                save_figure(g.figure, "{}_bg_remove".format(gene))

                gene = "IL7R"
                _umi = sc.get.obs_df(adata, keys = gene, layer='umi_counts')
                _bg = sc.get.obs_df(adata, keys = gene, layer='background_removed')
                data = pd.concat([_bg, _umi], axis=1)
                data.columns = ['cbender', 'umi']
                data = data.melt()
                g = sns.displot(data, x = 'value', hue = 'variable', kind = 'ecdf', log_scale=True)
                save_figure(g.figure, "{}_bg_remove".format(gene))

                # leiden on the tissue
                sc.tl.leiden(adata, key_added='{}_leiden'.format(tissue))
                sc.tl.rank_genes_groups(adata, groupby='{}_leiden'.format(tissue))
                sc.pl.rank_genes_groups(adata)
                
                # add meta data columns for B cell / not B cell
                
                ## Algo for definitely being a B cell
                ## TODO make leiden cluster cutoff depend to total number of clusters
                adata.obs["b_cell_super_cluster"] = adata.obs['my_leiden_majority_voting'].str.contains("B cells|Plasma")
                # change categorical to boolian
                bool_mapper = {"True":True, "False":False}
                adata.obs.predicted_doublets_umi_counts = adata.obs.predicted_doublets_umi_counts.astype(str).map(bool_mapper)
                adata.obs['probable_hq_single_b_cell'] = (adata.obs["b_cell_super_cluster"] & 
                                          ~adata.obs['predicted_doublets_umi_counts'] & 
                                          (~adata.obs['{}_leiden'.format(tissue)].astype(int) < 16))

                adata.obs["possible_b_cell"] = adata.obs["Immune_All_Low_predicted_labels"].str.contains("B cell|Plasma")

                adata.obs["probable_hq_single_not_b_cell"] = (~adata.obs["b_cell_super_cluster"] & 
                                              ~adata.obs['predicted_doublets_umi_counts']  
                                              & (~adata.obs['{}_leiden'.format(tissue)].astype(int) < 16))

                adata.obs["rare_or_bad_q_cell"] = adata.obs['{}_leiden'.format(tissue)].astype(int) > 16
                
                # Cell cycle Annotation
                cell_cycle_genes = pd.read_table('/home/michaelswift/repos/tabula-bursa/analysis/notebooks/cell_cycle_genes.tab', index_col=0)
                # develop a crude way to label cycling and non-cycling then test features by t-test
                same_length = False
                genes_contributing_to_score = 30
                if same_length:
                    sc.tl.score_genes(adata, gene_list=cell_cycle_genes.loc[:len(cycling), 'cc'], score_name='correlation_cycling')
                    sc.tl.score_genes(adata, gene_list=cell_cycle_genes.loc[:len(cycling), 'anti_cc'], score_name='anticorrelation_cycling')
                else:
                    sc.tl.score_genes(adata, gene_list=cell_cycle_genes.loc[:genes_contributing_to_score, 'cc'], score_name='correlation_cycling')
                    sc.tl.score_genes(adata, gene_list=cell_cycle_genes.loc[:genes_contributing_to_score, 'anti_cc'], score_name='anticorrelation_cycling')


                adata.obs['corr_cycling'] = adata.obs['correlation_cycling'] > 0.3
                adata.obs['anticorr_cycling'] = adata.obs['anticorrelation_cycling'] > 0.4

                cycling_mapper = {True:'True', False:'False'}

                use_non_cycling = False
                if use_non_cycling:
                    adata.obs['cycling'] = adata.obs['anticorr_cycling'].map(cycling_mapper)
                else:
                    adata.obs['cycling'] = adata.obs['corr_cycling'].map(cycling_mapper)
                
                fig, axs = plt.subplots(1,2, sharey=True)
                fig.tight_layout(pad=2.0)
                x, y = adata.obs['anticorrelation_cycling'], adata.obs['correlation_cycling']
                sns.histplot(x, ax=axs[0])
                sns.histplot(y, ax=axs[1])
                plt.vlines(x = 0.1, ymax=10000, ymin = 0, linestyles='dotted', color = 'r')
                plt.yscale('log')
                save_figure(fig, "{}_cell_cycle_scoring".format(tissue))
                # write obs matrix for vdj integration
                adata.write_h5ad("outputs/tissue_objs/{}_annotated_processed.h5ad".format(tissue))
                df = adata.obs[['sample_uid', 'Immune_All_Low_predicted_labels', 'Immune_All_High_predicted_labels', 
                                'Immune_All_High_conf_score', 'Immune_All_Low_conf_score', 'doublet_scores_umi_counts', 'predicted_doublets_umi_counts', 
                                'n_genes', 'total_counts', 'total_counts_mt', '{}_leiden'.format(tissue), 'tissue', 'probable_hq_single_b_cell', "possible_b_cell", "probable_hq_single_not_b_cell", "rare_or_bad_q_cell"]]
                df.to_csv('{}_gex_labels.tab'.format(tissue), sep = '\t')
                # perform garbage collection:
                del(adata)
                gc.collect()



# run routine
tissues = ["LN", "BM", "SP", "PB"]
for h5ad_name in snakemake.input
h5ad_name = "/home/michaelswift/repos/shared_data/pipeline_outs/"

run_qc_routine(data_path, tissues)

dfs = []
for tissue in tissues:
    df = pd.read_table('{}_gex_labels.tab'.format(tissue), index_col=0)
    dfs.append(df)
df = pd.concat(dfs)
df.to_csv(gex_labels, sep = '\t')
