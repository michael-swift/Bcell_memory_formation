### Dictionaries for color definition etc.
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from matplotlib import pyplot as plt
import seaborn as sns
import scanpy as sc
import scanpy.external as sce
from matplotlib.lines import Line2D

############################
##### Scanpy Helper Fxns ###
############################

def cluster(adata, batch_correct=False, batch_key="tissue", scale = False, pca = True):
    if pca:
        print("PCA-ing")
        sc.pp.pca(adata)
    print("drawing neighbor graph")
    if batch_correct == True:
        print("batch corrected neighbors")
        sce.pp.bbknn(adata, batch_key=batch_key)
    else:
        sc.pp.neighbors(adata, n_neighbors=20)
    print("UMAP-ing")
    sc.tl.umap(adata)
    print("leiden-ing")
    sc.tl.leiden(adata, resolution=0.2)
    return adata

def umap_color_1(adata, column_name, variable_name, **kwargs):
    adata.obs['{}_cluster'.format(variable_name)] = (adata.obs[column_name] == variable_name).astype(str)
    sc.pl.umap(adata, color = "{}_cluster".format(variable_name), **kwargs)
    return

def reformat_adata_obs(adata):
    df = adata.obs.astype('object')

    cols = df.columns
    for c in cols:
        try:
            df[c] = pd.to_numeric(df[c])
        except:
            pass
    adata.obs = df
    return adata


def safe_subset(adata, subsetting_var, status, group):
    subset_adata = adata[adata.obs[subsetting_var] == status]
    subset_adata = subset_adata[subset_adata.obs[group].dropna().index]
    # filter out categories where there aren't enough cells to perform DE
    selector = subset_adata.obs[group].value_counts() < 10
    selector = list(selector[selector == False].index)
    subset_adata = subset_adata[subset_adata.obs[group].isin(selector)]
    return subset_adata

def DE_analysis(adata, subsetting_var, status, group, min_logfoldchange=1):
    subset_adata = safe_subset(adata, subsetting_var, status, group)
    #subset_adata = adata
    key_tag = status + "_" + group

    sc.tl.rank_genes_groups(
        subset_adata, groupby=group, method="wilcoxon", key_added=key_tag
    )
    sc.pl.rank_genes_groups(subset_adata, key=key_tag, save="_{}".format(key_tag))
    sc.tl.dendrogram(subset_adata, groupby=group)
    sc.pl.rank_genes_groups_stacked_violin(
        subset_adata,
        key=key_tag,
        n_genes=10, min_logfoldchange=min_logfoldchange,
        save="_{}".format(key_tag),
    )
    sc.pl.rank_genes_groups_matrixplot(
        subset_adata,
        key=key_tag,
        n_genes=12, min_logfoldchange=min_logfoldchange,
        save="_{}".format(key_tag),
    )
    
    sc.pl.rank_genes_groups_matrixplot(
        subset_adata,
        key=key_tag,
        n_genes=12, min_logfoldchange=min_logfoldchange,
        save="_{}".format(key_tag), values_to_plot = 'logfoldchanges'
    )
    sc.pl.rank_genes_groups_matrixplot(
        subset_adata,
        key=key_tag,
        n_genes=12, min_logfoldchange=min_logfoldchange,
        save="_{}".format(key_tag), values_to_plot = 'pvals_adj'
    )
    result = subset_adata.uns[key_tag]
    groups = result["names"].dtype.names
    # construct dataframe of top 100 DE genes
    df = pd.DataFrame(
        {
            group + "_" + key[:1]: result[key][group]
            for group in groups
            for key in ["names", "pvals_adj"]
        }
    ).head(300)
    # Use top genes to create a group score
    n_top_genes = 20
    score_names = []
    print(subset_adata.obs[group].unique())
    for cat in subset_adata.obs[group].unique():
        print(cat)
        construct = cat + "_n"
        score_name = cat + "_score"
        sc.tl.score_genes(
            subset_adata,
            gene_list=df.loc[:n_top_genes, construct],
            score_name=score_name,
        )
        score_names.append(score_name)

    for color in score_names + UMAP_variables_of_interest:
        sc.pl.umap(subset_adata, color=color, save="_{}_{}".format(color, key_tag))
    return subset_adata, df

###########################
## Plotting Dictionaries ##
###########################
def celltype_colors():
    color_dict = {'Plasma B cell': "#29975D" , 'Myeloid cell': "#972963",
                    'Resting B cell': "#6BB9CB",
       'Activated B cell': "#000000", 'T cell': "#EF9A9A", 'NK / Dead Cell':"#ABABAB" , 'NK cell': "#8A679A"}
    return color_dict


def tissue_colors():
    tissue_color_dict = {"PB": "#d6616b", "BM":"#cedb9c", "LN":"#8c6d31", "SP":"#393b79"}
    return tissue_color_dict

def sample_id_colors():
    sample_id_color_dict = {'BM CD138+': '#1f09ff',
             'BMMNC': '#0977ff',
             'Day 0': '#14e0be',
                            'Day 12':"#ff3939",
                            'Day 8':"#da5435",
                            'Day 4':"#c0cd30",
                            'PBMC':"#09a3ff"                         }
    return sample_id_color_dict

def switched_colors():
    color_dict = {"IGHM|IGHD": "#14e0be",
                  "switched":"#000000"}
    return color_dict

def mutation_colors():
    mutation_colors = {'mutated': '#FF5733', 'germline' : '#004D40', 'heavily mutated':'#581845'}
    return mutation_colors

def IGH_colors():
    IGH_color_dict = {'IGHA1': '#FFAB91',
             'IGHA2': '#FFD180',
             'IGHD': '#00E676',
             'IGHM': '#69F0AE',#"#00ff7f",
             'IGHG2': '#B388FF',
             'IGHG4': '#008784',
             'IGHG1': '#0B3954',
             'IGHG4': '#6A1B9A',
             'IGHG3': '#536DFE',
             'IGHE':'#FF5A5F'}
    return IGH_color_dict

def IGH_simplify():
    IGH_simple = {'IGHA1': 'IGHA',
             'IGHA2': 'IGHA',
             'IGHD': 'IGHD',
             'IGHM': 'IGHM',#"#00ff7f",
             'IGHG2': 'IGHG',
             'IGHG4': 'IGHG',
             'IGHG1': 'IGHG',
             'IGHG4': 'IGHG',
             'IGHE':'IGHE'}
    return IGH_simple

def IGH_colors_simple():
    
    IGH_color_dict = {'IGHA': '#e0b96a',
             'IGHD': '#798658',
             'IGHM': '#17a589',#"#00ff7f",
             'IGHG': '#0D47A1',
             'IGHE':'#8b4cf4'}
    return IGH_color_dict

def plot_order():
    order = ['BMMNC',
 'PBMC',
 'BM CD138+',
 'Day 0',
 'Day 4',
 'Day 8',
 'Day 12']
    return order

def timecourse():
    order = [
 'Day 0',
 'Day 4',
 'Day 8',
 'Day 12']
    return order

def pseudo_timecourse():
    order = [
 'Day 0',
 'Day 4',
 'Day 8',
 'Day 12', 'BM CD138+']
    return order

mutation_cutoffs = [0.02, 0.10]

def IGH_switched():
    
    IGH_switched = {'IGHA1': 'switched',
             'IGHA2': 'switched',
             'IGHD': 'IGHM|D',
             'IGHM': 'IGHM|D',#"#00ff7f",
             'IGHG2': 'switched',
             'IGHG4': 'switched',
             'IGHG1': 'switched',
             'IGHG4': 'switched',
             'IGHE':'switched'}
    return IGH_switched

def timepoint_shapes():
    timepoint_shapes_dict = {"Day 0": "triangle", "Day 4": "square", "Day 8": "pentagon", "Day 12": "octagon"}
    return timepoint_shapes_dict

def celltype_colors():
    color_dict = {'Plasma B cell': "#29975D" , 'Myeloid cell': "#972963",
                    'Resting B cell': "#6BB9CB",
       'Activated B cell': "#000000", 'T cell': "#EF9A9A", 'NK / Dead cell':"#ABABAB" , 'NK cell': "#8A679A"}
    return color_dict


def bcelltype_colors_dict():
    colors = {'B cells' : '#8000ff', 'Memory B cells' : '#1996f3', 'Naive B cells':'#4df3ce', 'Plasma cells':'#b2f396', 'Plasmablasts': '#ff964f', 'Prolif. GC B cells': '#ff0000'}
    return colors


IGH_order = {'IGHM':0, 'IGHD':1, 'IGHG3':2, 'IGHG1':3, 'IGHA1':4, "IGHG2":5, 'IGHG4':6, "IGHE":7, "IGHA2":8}

IGH_simple_markers = {'IGHM':'^', 'IGHD':"v", 'IGHG':"o", 'IGHA':"s","IGHE":"*"}
UMAP_variables_of_interest = ['sample_id', 'isotype_simple', 'switched', 'simple_mutation_status']


###########################
##        Execute        ##
###########################
# loads dictionaries into memory
timepoint_shapes = timepoint_shapes()
celltype_colors = celltype_colors()
bcelltype_colors = bcelltype_colors_dict()
sample_id_colors = sample_id_colors()
switched_colors = switched_colors()
mutation_colors = mutation_colors()
plot_order = plot_order()
timecourse = timecourse()
pseudo_timecourse = pseudo_timecourse()
igh_colors = IGH_colors()
igh_colors_simple = IGH_colors_simple()
igh_genes = list(IGH_order.keys())
# other settings
savefig_args = {"dpi": 300, "bbox_inches": "tight", "pad_inches": 0, "transparent": True}
output_suffix = ""
output_formats = [".png", ".svg"]


###########################
##  Seqclone Helper Fxns ##
###########################
def selection_helper(df, param, value):
    selector = df[param].value_counts() > value
    idxs = selector[selector == True].index
    return(df[df[param].isin(idxs)])

def shuffleWithinGroup(data, label, group_name):
    """ returns a shuffled series of the label (column name), shuffled within the groups supplied """
    list_of_dfs = []
    for group, frame in data.groupby(group_name):
        frame.loc[:,label] = np.random.permutation(frame.loc[:,label])
        list_of_dfs.append(frame)
    df = pd.concat(list_of_dfs)
    return df.loc[:,label]

# visualize the color dictionaries
def plot_colortable(colors, title, sort_colors=True, emptycols=0):

    cell_width = 230
    cell_height = 20
    swatch_width = 48
    margin = 20
    topmargin = 40

    # Sort colors by hue, saturation, value and name.
    by_hsv = [(v, k) for k, v in colors.items()]
    
    if sort_colors is True:
        by_hsv = sorted(by_hsv)
    names = [name for hsv, name in by_hsv]

    n = len(names)
    ncols = 2 - emptycols
    nrows = n // ncols + int(n % ncols > 0)

    width = cell_width * 4 + 2 * margin
    height = cell_height * nrows + margin + topmargin
    dpi = 72

    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin/width, margin/height,
                        (width-margin)/width, (height-topmargin)/height)
    ax.set_xlim(0, cell_width * 4)
    ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()
    ax.set_title(title, fontsize=24, loc="left", pad=10)

    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height

        swatch_start_x = cell_width * col
        swatch_end_x = cell_width * col + swatch_width
        text_pos_x = cell_width * col + swatch_width + 7

        ax.text(text_pos_x, y, name, fontsize=12,
                horizontalalignment='left',
                verticalalignment='center')

        ax.hlines(y, swatch_start_x, swatch_end_x,
                  color=colors[name], linewidth=18)

    return fig
