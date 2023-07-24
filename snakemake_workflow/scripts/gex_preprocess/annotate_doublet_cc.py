import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import gc
import pathlib
import seaborn as sns
from matplotlib import pyplot as plt
import anndata as ad

###############################
###############################
# single cell preprocessing ###
######  functions #############
###############################
###############################
def cluster(adata):
    sc.pp.neighbors(adata, use_rep="X_scVI_all")
    sc.tl.leiden(adata, resolution = 1.5)
    return adata

def assign_technical_categories(adata):
    """Possible Categories are
    probable_hq_single_b_cell: Highest Confidence Annotation, low risk of FP
    probable_hq_single_not_b_cell: Highest Confidence Annotation, low risk of FP
    possible_b_cell: This is likely a high quality cell or doublet, it is likely to be a B cell
    rare_or_bad_q_cell: empirically these seem to be non-B cells and likely ambient droplets or other technical artifcats
    """
    ## Flag Believable leidens
    believable_leidens = adata.obs['leiden'].value_counts(normalize=True)
    believable_leidens = believable_leidens[believable_leidens > 0.001].index.tolist()
    ## doublet flagging by majority voting of high hierarchy cell type
    majority_voting_doublet = pd.Series(False, index=adata.obs.index)
    # iterate over each unique leiden value
    for i in adata.obs.leiden.unique():
        max_value = adata.obs.groupby('leiden')['Immune_All_High_predicted_labels'].value_counts(normalize = True).xs(i).max()
        if max_value < 0.90:
            # assign True to cells in this leiden group
            majority_voting_doublet[adata.obs['leiden'] == i] = True
    # add the majority_voting_doublet column to adata.obs
    adata.obs['majority_voting_doublet'] = majority_voting_doublet
    adata.obs["rare_or_bad_q_cell"] = ~adata.obs["leiden"].isin(believable_leidens)
    # manually change majority voting flag for plasmablasts
    condition = adata.obs["Immune_All_Low_predicted_labels"] == "Plasmablasts"
    # Updating the "majority_voting_doublet" column based on condition
    adata.obs.loc[condition, "majority_voting_doublet"] = False

    ##########################################
    # Rules to Assign Technical Categories
    adata.obs["possible_b_cell"] = (
        adata.obs["majority_voting_low_predicted_labels"].str.contains("B cell|Plasma")
        | adata.obs["Immune_All_High_predicted_labels"].str.contains("B cell|Plasma")
        | adata.obs["Immune_All_Low_predicted_labels"].str.contains("B cell|Plasma")
    )
    # change categorical doublet label to boolean for filtering
    bool_mapper = {"True": True, "False": False}
    adata.obs["predicted_doublets_counts"] = adata.obs["predicted_doublets_counts"].astype(str).map(bool_mapper)
    adata.obs["probable_hq_single_b_cell"] = (
        adata.obs["possible_b_cell"]
        & (adata.obs["Immune_All_Low_predicted_labels"].str.contains("B cell|Plasma"))
        & (~adata.obs["majority_voting_doublet"])
           & (adata.obs['Immune_All_Low_conf_score'] > 0.95))
    adata.obs["probable_hq_single_not_b_cell"] = (
        (~adata.obs["possible_b_cell"])
        & (~adata.obs["rare_or_bad_q_cell"])
        & (~adata.obs["majority_voting_doublet"])
        & (adata.obs['Immune_All_Low_conf_score'] > 0.95)
    )
    return adata

def annotate_cell_cycle(adata, cell_cycle_genes, genes):
    genes_contributing_to_score = genes
    sc.tl.score_genes(
        adata,
        gene_list=cell_cycle_genes.loc[:genes_contributing_to_score, "cc"],
        score_name="correlation_cycling",
    )
    sc.tl.score_genes(
        adata,
        gene_list=cell_cycle_genes.loc[:genes_contributing_to_score, "anti_cc"],
        score_name="anticorrelation_cycling",
    )
    adata.obs["corr_cycling"] = adata.obs["correlation_cycling"] > 0.01
    adata.obs["anticorr_cycling"] = adata.obs["anticorrelation_cycling"] > 0.4
    # Make Categories non-boolean for plotting in Scanpy
    cycling_mapper = {True: "True", False: "False"}
    use_non_cycling = False
    if use_non_cycling:
        adata.obs["cycling"] = adata.obs["anticorr_cycling"].map(cycling_mapper)
    else:
        adata.obs["cycling"] = adata.obs["corr_cycling"].map(cycling_mapper)
    return adata

def convert_non_numerical_and_boolean_to_str(df):
    for col in df.columns:
        if pd.api.types.is_categorical_dtype(df[col].dtype):
            df[col] = df[col].astype(str)
        elif not np.issubdtype(df[col].dtype, np.number) or np.issubdtype(df[col].dtype, np.bool_):
            df[col] = df[col].astype(str)
    return df

###############################
#### execution ################
###############################

h5ad = str(snakemake.input)
output_file = str(snakemake.output)
cell_cycle_genes = pd.read_table(str(snakemake.params.cell_cycle_genes), index_col=0)
adata = sc.read_h5ad(h5ad)
# run routine
adata = cluster(adata)

    # Cell Cycle Annotation
adata = annotate_cell_cycle(adata, cell_cycle_genes, genes = 30)

adata = assign_technical_categories(adata)

# convert to strings which should avoid h5py error?
adata.obs = convert_non_numerical_and_boolean_to_str(adata.obs)

print("writing h5ad")
adata.write_h5ad(
    output_file,
    compression="gzip",

)
