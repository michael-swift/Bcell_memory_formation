import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt 

import time

#tree stuff
from Bio import Phylo
from io import StringIO
import pylab
from matplotlib.patches import Circle, PathPatch


############## FUNCTIONS ######################

def get_x_positions(tree):
        """Create a mapping of each clade to its horizontal position.
        Dict of {clade: x-coord}
        """
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        return depths

def get_y_positions(tree):
        """Create a mapping of each clade to its vertical position.
        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        maxheight = tree.count_terminals()
        # Rows are defined by the tips
        heights = dict((tip, maxheight - i)
                       for i, tip in enumerate(reversed(tree.get_terminals())))

        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (heights[clade.clades[0]] +
                              heights[clade.clades[-1]]) / 2.0

        if tree.root.clades:
            calc_row(tree.root)
        return heights


def get_label(leaf):
        return leaf.name


def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
                names[clade.name] = clade
    return names
    
    
def draw_tree(target_lineage, df, tree_df, color_by='tissue', color_dict=None):
    
    # set up tree data structure
    lineage_df = df[df.lineage_uid == target_lineage]
    tree_nwk = tree_df.loc[target_lineage, 'v_phylogeny']


    vseq_df=lineage_df[['v_seq_id','v_mismatch',color_by,'v_seq_len','total_reads','total_umis']]#.drop_duplicates(ignore_index=True)

    vseq_df = vseq_df.groupby(['v_seq_id',
                               'v_mismatch',
                               'v_seq_len'])[[color_by,
                                            'total_reads',
                                            'total_umis']].agg({color_by:set,
                                                                'total_reads':sum,
                                                                'total_umis':sum}).reset_index()
    vseq_df = vseq_df.rename(columns={0:color_by})
    multi_tissue = vseq_df[color_by].str.len() > 1

    vseq_df.loc[multi_tissue,color_by] = 'multiple'
    vseq_df[color_by] = vseq_df[color_by].apply(lambda x: x if x == 'multiple' else list(x)[0])

    vseq_df=vseq_df.reset_index().set_index('v_seq_id')

    subsets = sorted(vseq_df[color_by].unique())

    tree = Phylo.read(StringIO(tree_nwk), 'newick')
    
    

    it = 0

    #root on clade closest to germline

    min_d=300
    idx_min=-1
    closest_clade = None
    for clade in tree.find_clades():
        if clade.name:
            if not('inferred' in clade.name):
                d = vseq_df.loc[clade.name,'v_mismatch'].astype(int)
            else:
                d = 0
            if d < min_d:
                min_d = d
                idx_min = clade.name
                closest_clade = clade


    if min_d == 0:
        germ_clade = closest_clade
    else:
        germ_clade = None

    tree.root_with_outgroup(closest_clade)
    tree.ladderize()
    tree.rooted = True
    
    # set tree style
    tree_color = 'k'
    tree.root.color= tree_color
    tree.root.width = 0.3
    fig_height = 3*len(tree.get_terminals())/100
    fig_width = 10

    if fig_height>10:
        fig_height=10
    fig, ax = pylab.subplots(figsize=(fig_width,fig_height))

    
    x_positions = get_x_positions(tree)
    y_positions = get_y_positions(tree)

    if color_dict is None:
        subset_colors = {}
    else:
        subset_colors = color_dict
    
    for feature in subsets:
        if feature in subset_colors.keys():
            line1, = ax.plot([],[],lw=8,label=feature, color=subset_colors[feature])
        else:
            line1, = ax.plot([],[],lw=8,label=feature)
            subset_colors.update({feature:line1.get_color()})
    subset_colors['multiple'] = 'k'

    for clade_id in x_positions.keys():
        if (clade_id.name is None) or (clade_id.name.startswith("Inner")):
            pass
        else:
            x = x_positions[clade_id]
            y = y_positions[clade_id]
            if clade_id.name.endswith("inferred"):
                kwargs = dict(color='0.5',
                              lw=.5,
                              edgecolor=(0,0,0,1.),
                              zorder=10,
                              s=12,
                              alpha=.6)
            elif (clade_id.name == idx_min) and (vseq_df.loc[clade_id.name, 'v_mismatch'] == 0):
                clade_subset = vseq_df.loc[clade_id.name, color_by]
                s = 6+12* ((np.log(vseq_df.loc[clade_id.name, 'total_umis'])/np.log(10))**2)
                error_corrected = (vseq_df.loc[clade_id.name, 'total_reads'] > 
                                    vseq_df.loc[clade_id.name, 'total_umis']) and (vseq_df.loc[clade_id.name, 'total_reads'] > 2)
                if vseq_df.loc[clade_id.name, 'total_reads'] > vseq_df.loc[clade_id.name, 'total_umis']:
                    marker = 'o'
                    alpha = 0.6
                else:
                    marker = 'o'
                    alpha = 0.6

                s = 12
                kwargs = dict(color=subset_colors[clade_subset],
                               marker=marker,
                              lw=.5,
                              edgecolor=(0,0,0,1.),
                              zorder=10,
                              s=s,
                              alpha=alpha)
            else:
                clade_subset = vseq_df.loc[clade_id.name,
                          color_by]
                s = 6+12* ((np.log(vseq_df.loc[clade_id.name, 'total_umis'])/np.log(10))**2)
                if vseq_df.loc[clade_id.name, 'total_reads'] > vseq_df.loc[clade_id.name, 'total_umis']:
                    marker = 'o'
                    alpha = 0.6
                else:
                    marker = 'o'
                    alpha = 0.6

                s = 12
                kwargs = dict(color=subset_colors[clade_subset],
                              lw=0.,
                              marker=marker,
                              edgecolor=tree_color,
                              zorder=10,
                              s=s,
                              alpha=alpha)

            ax.scatter(x,y, clip_on=False, **kwargs)


    Phylo.draw(tree,
           label_func=get_label,
           do_show=False,
           show_confidence=False, 
           axes = ax
           )
    ax.axis('off')
    ax.set_ylim([-16,ax.get_ylim()[0]])

    ax.set_xlim([-0.05,0.5])
    ax.plot([0.00, 0.05], [-10,-10], color = 'k')
    ax.legend(loc='upper left', frameon=False)
    
    return fig, ax
