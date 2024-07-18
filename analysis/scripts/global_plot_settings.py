import matplotlib as mpl
import numpy as np
import pandas as pd

# Define hex codes and color palettes
hex_colors = {
    'dark_purple': "#6F4B70",
    'grey_green': "#abb8b8",
    'sky_blue': "#78BACC",
    'green1': "#77BA99",
    'green2': "#68918f",
    'green3': "#266967",
    'thistle': "#DAC1E0"
}

# Convert hex to RGB
rgb_colors = [mpl.colors.to_rgb(color) for color in hex_colors.values()]
cdict = {'red': [x[0] for x in rgb_colors],
         'green': [x[1] for x in rgb_colors],
         'blue': [x[2] for x in rgb_colors]}

# Create custom colormap
donor_cmap = mpl.colors.LinearSegmentedColormap.from_list('donor_cmap', list(hex_colors.values()), N=6)
donor_cycle = [donor_cmap(i/6) for i in range(6)]

# Color dictionaries
donor_colors = {f'TBd{i+1}': color for i, color in enumerate(donor_cycle)}
tissue_colors = {'SP': '#801630FF', 'PB': '#C1717B', 'BM': '#E05353', 'LN': '#E4B363', 'both': 'k', 'multiple': 'k'}
IGH_colors = {'IGHA1': '#fb9a99', 'IGHA2': '#e31a1c', 'IGHD': '#fec44f', 'IGHM': '#33a02c', 'IGHG2': '#02818a', 'IGHG4': '#6a51a3', 'IGHG1': '#386cb0', 'IGHG3': '#3690c0', 'IGHE': '#67001f', "ambiguous": '#f7f7f7'}
bcelltype_colors = {'Memory B cells': '#5cecffff', 'Naive B cells': '#ff61c6ff', 
                    'Plasma cells': '#b2f396', 'Plasmablasts': '#ff9900ff', 
                    'Proliferative Germinal Centre B cells': '#0a0c37ff', 'Age-associated B cells': '#375971ff'}

bcelltype_colors_alt = {
                    'Plasmablasts':'#a6cee3',
                    'Plasma cells':'#1f78b4',
                    'Naive B cells':'#b2df8a',
                    'GC B cells':'#33a02c',
                    'ABCs':'#fb9a99',
                    'Memory B cells':'#e31a1c',
                    'ASC-1':'#A6CEE3',
                    'ASC-2':'#4598C4',
                    'ASC-3':'#1B699D',
                    'ASC-4':'#7071A3',
                    'ASCs': '#1f78b4',
                    'B cells': "#e31a1c",
                    'Pro-B cells':'#FFFFFF',
                    'Cycling B cells':'#FFFFFF'}

asc_colors = {
    'ASC-1': '#A6CEE3',
    'ASC-2': '#4598C4',
    'ASC-3': '#1B699D',
    'ASC-4': '#7071A3'
}

memory_b_cell_colors = {
    'NS-Memory': '#fdb863',
    'NS-ABCs': '#e08214',
    'NS-Memory-CD21++': '#b35806',
    'SW-Memory': '#b2abd2',
    'SW-ABCs': '#8073ac',
    'SW-Memory-CD21++': '#542788'
}


# Simplified mappings
celltypist_simpler = {
    "Proliferative germinal center B cells": "GC B cells",
    "Germinal center B cells": "GC B cells",
    "Age-associated B cells": "ABCs",
    "Plasmablasts": "ASCs",
    "Plasma cells": "ASCs"
}

def IGH_switched():
    
    IGH_switched = {'IGHA1': 'SW',
             'IGHA2': 'SW',
             'IGHD': 'NS',
             'IGHM': 'NS',
             'IGHG2': 'SW',
             'IGHG3': 'SW',
             'IGHG1': 'SW',
             'IGHG4': 'SW',
             'IGHE':'SW'}
    return IGH_switched

# Helper functions
def get_colors(color_type):
    color_dicts = {
        'donor': donor_colors,
        'tissue': tissue_colors,
        'IGH': IGH_colors,
        'bcelltype': bcelltype_colors,
        'bcelltype_alt': bcelltype_colors_alt,
        'asc': asc_colors,
        'memory_b_cell': memory_b_cell_colors
    }
    return color_dicts.get(color_type, {})

def set_colors(color_type, new_dict):
    if color_type in ['donor', 'tissue', 'IGH', 'bcelltype']:
        globals()[f"{color_type}_colors"].update(new_dict)

def convert_boolean_columns(df):
    for col in df.columns:
        unique_values = df[col].dropna().unique()
        if set(unique_values) <= {'True', 'False', 'nan'}:
            df[col] = df[col].replace({'True': True, 'False': False, 'nan': np.nan}).astype('boolean')
    return df