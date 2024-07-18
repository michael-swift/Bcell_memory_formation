import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import matplotlib.patches as patches
import packcircles as pc


def assign_colors(keys, 
                  cmap=None, 
                  vmax=None, 
                  vmin=None):
    n_colors = len(keys)

    if not pd.api.types.is_numeric_dtype(pd.Series(keys)):     
        if cmap is None:
            cmap='Set2'
        x = np.arange(0,1,n_colors)
            
    else:
        if cmap is None:
            cmap='mako'
        if vmin is None:
            vmin=min(keys)
        if vmax is None:
            vmax=max(keys)
        norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        x = norm(keys)

    colors = mpl.colormaps[cmap](x)#.to_rgba(x)
    return {keys[i]:colors[i] for i in range(n_colors)}


def pack_circles(radii, col_names=['x','y','r'], sortx=False):
    circles = pc.pack(radii)
    df = pd.DataFrame(list(circles))
    df.columns = col_names
    if sortx:
        df = df.sort_values([col_names[2], col_names[0],col_names[1]], ignore_index=True)
    else:
        df = df.sort_values([col_names[2]], ignore_index=True, ascending=False)
    return df

def get_within_group_coordinates(n, r=2, pad=3):
    if n > 2:
        radii = np.ones(n)*r
        df = pack_circles(radii, col_names=['dx','dy','r'], sortx=True)
        max_dist = np.max(np.sqrt(df.dx**2 + df.dy**2))
        R = max_dist + r + pad
    elif n==1:
        df = pd.DataFrame({'dx':[0], 'dy':[0], 'r':[r]})
        R = r + pad
    elif n==2:
        df = pd.DataFrame({'dx':[-r,r], 'dy':[0,0], 'r':[r,r]})
        R = 2*r + pad
    return df, R

def get_hex_coordinates(data_frame,
                        group_col='lineage_id',
                        sort=None,
                        hue=None,
                        palette={},
                        r=2, 
                        pad=3,
                        max_colors=10,
                        cmap=None
                        ):
    
    """ Computes packing coordinates within group """
    group_df = data_frame.groupby(group_col).size().reset_index().rename(columns={0:'n'})
    unique_group_sizes = group_df['n'].unique()

    # first compute packing coordinates within groups
    ## this will make it possible to compute group radii,
    ## which can then be use to compute group centroid coordinates
    within_group_dict = {n: get_within_group_coordinates(n,r)
                        for n in unique_group_sizes
                        }

    # compute centroid coordinates
    group_df['R_c'] = group_df['n'].map(lambda x: within_group_dict[x][1])
    group_df = group_df.sort_values('R_c', ignore_index=True, ascending=False)

    group_centroids = pack_circles(group_df['R_c'].values, col_names=['x_c','y_c','R_c'])
    group_df = group_df.merge(group_centroids, left_index=True, 
                                               right_index=True)
    if (group_df['R_c_x'] == group_df['R_c_y']).all():
        group_df = group_df.rename(columns={'R_c_x':'R_c'})
        group_df = group_df.drop('R_c_y', axis=1)
    group_df = group_df.set_index(group_col)

    # assign coordinates to individual objects
    
    coord_df = data_frame.copy()
    
    coord_df['x_c'] = coord_df[group_col].map(group_df['x_c'])
    coord_df['y_c'] = coord_df[group_col].map(group_df['y_c'])
    coord_df['R_c'] = coord_df[group_col].map(group_df['R_c'])
    coord_df['n_group'] = coord_df[group_col].map(group_df['n'])
    
    if sort is None:
        sort_key=hue
    else:
        sort_key=sort
    # add offsets from centroid
    if hue is None:
        sort_keys = [group_col]
    else:
        sort_keys = [group_col, sort_key]
        
    coord_df = coord_df.sort_values(sort_keys, ignore_index=True)

    coords = ['x', 'y']
    coord_df['r'] = r
    coord_df[coords] = np.nan
    for group_id in coord_df[group_col].unique():
        subset = coord_df[group_col] == group_id
        n = group_df.loc[group_id,'n']

        for coord in coords:
            coord_df.loc[subset, coord] = coord_df[f'{coord}_c'] 
            coord_df.loc[subset, coord] += within_group_dict[n][0][f'd{coord}'].values

    # finally, assign colors
    if hue is None:
        #color by group
        hue=group_col
        group_keys = {group_df.index[i]: i%max_colors for i in np.arange(group_df.shape[0])}
        group_colors = assign_colors(np.arange(min(max_colors,group_df.shape[0])), cmap=cmap)
        group_colors = {k:group_colors[v] for k,v in group_keys.items()}
        
    else:
        group_colors=palette
        #  check if palette contains all colors
        unique_keys = data_frame[hue].unique()
        keys_present = [key in palette.keys() for key in unique_keys]
        if all(keys_present):
            pass
        else:
            missing_keys = unique_keys[~np.asarray(keys_present, dtype=bool)]
            additional_colors = assign_colors(missing_keys, cmap=cmap)
            group_colors.update(additional_colors)
            
    coord_df['color'] = coord_df[hue].map(group_colors)
    
    return coord_df

def plot_hexplot(data_frame, 
                 fig_ax=None,
                 group_col='lineage_id',
                 sort=None,
                 hue=None,
                 palette={},
                 r=2, 
                 pad=1,
                 max_colors=10,
                 cmap=None):
    
    if fig_ax is None:
        fig, ax = plt.subplots()
        
    
    data_frame = get_hex_coordinates(data_frame, 
                                     group_col=group_col, 
                                     hue=hue,
                                     sort=sort,
                                     palette=palette,
                                     r=r,
                                     pad=pad,
                                     max_colors=max_colors,
                                     cmap=cmap)

    y_max, y_min = (data_frame.y +  data_frame.r).max(), (data_frame.y -  data_frame.r).min()
    x_max, x_min = (data_frame.x +  data_frame.r).max(), (data_frame.x -  data_frame.r).min()
    y_range = y_max - y_min
    x_range = x_max - x_min

    fig_height = y_range/100*2
    fig_width = x_range/100*2
    fig.set_figheight(fig_height)
    fig.set_figwidth(fig_width)

    for idx, row in data_frame.iterrows():
        patch = patches.CirclePolygon(
            (row['x'],row['y']),
            row['r'], resolution=6,
            color=row['color'],
            alpha=1
        )
        ax.add_patch(patch)

    ax.set(xlim=(x_min, x_max), ylim=(y_min, y_max))

    plt.axis('off')
    return fig, ax, data_frame
