import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
import matplotlib as mpl

class ColorDictionaries:
    def __init__(self):
        self.color_dicts = {
            "switched_colors": {"IGHM|IGHD": "#14e0be", "switched":"#000000"},
            "mutation_colors": {'mutated': '#FF5733', 'germline' : '#004D40', 'heavily mutated':'#581845'},
            "IGH_colors": {
                'IGHA1': '#fb9a99',
                'IGHA2': '#e31a1c',
                'IGHD': '#fec44f',
                'IGHM': '#33a02c',
                'IGHG2': '#02818a',
                'IGHG4': '#6a51a3',
                'IGHG1': '#386cb0',
                'IGHG3': '#3690c0',
                'IGHE':'#67001f',
                "ambiguous":'k'
            }, "IGH_simple_colors": {'IGHA': '#e31a1c',
             'IGHD': '#fec44f',
             'IGHM': '#33a02c',
             'IGHG': '#386cb0',
             'IGHE':'#67001f',
                "ambiguous":'#f7f7f7'     }, 
        "all_bcelltype_colors": {'Memory B cells' : '#5cecffff', 'Naive B cells' : '#ff61c6ff', 'Plasma cells':'#b2f396',
               'Plasmablasts': '#ff9900ff', 'Proliferative Germinal Centre B cells': '#0a0c37ff',
                 'Age-associated B cells': '#375971ff', 'Transitional B cells': "#b3b5bb", 
                 'Large pre-B cells': '#b3b5bb', 'Pro-B cells': '#b3b5bb', 'Small pre-B cells': 
                 '#b3b5bb', 'Pre-pro-B cells':'#b3b5bb', 'Proliferative germinal center B cells':'#0a0c37ff', 
                 'B cells':'#b3b5bb', 'Germinal center B cells': '#817f75', 'Cycling B cells':'#b3b5bb'},
        "bcelltype_colors_full" :{
                    'Plasmablasts':'#a6cee3',
                    'Plasma cells':'#1f78b4',
                    'Naive B cells':'#b2df8a',
                    'Germinal center B cells':'#fdbb84',
                    'Age-associated B cells':'#fb9a99',
                    'Memory B cells':'#e31a1c',
                    'Plasmablasts':'#a6cee3',
                    'Plasma cells':'#1f78b4',
                    'Naive B cells':'#b2df8a',
                    'Proliferative germinal center B cells':'#33a02c',
                    'Age-associated B cells':'#fb9a99',
                    'Memory B cells':'#e31a1c'},
       "bcelltype_colors_alt":{
                    'Plasmablasts':'#a6cee3',
                    'Plasma cells':'#1f78b4',
                    'Naive B cells':'#b2df8a',
                    'GC B cells':'#33a02c',
                    'ABCs':'#fb9a99',
                    'Memory B cells':'#e31a1c'},
        "asc_colors": {'ASC-1':'#A6CEE3',
                                         'ASC-2':'#4598C4',
                                         'ASC-3':'#1B699D',
                                         'ASC-4':'#7071A3'},
        "memory_b_cell_colors": {
                    'NS-Memory':'#fdb863',
                     'NS-Memory-ABC':'#e08214',
                      'NS-Memory-CD21++':'#b35806',
                    'SW-Memory':'#b2abd2',
                    'SW-Memory-ABC':'#8073ac',
                    'SW-Memory-CD21++':'#542788'}}

    def get_colors(self, dict_name):
        return self.color_dicts.get(dict_name, {})

celltypist_simpler = {
    "Proliferative germinal center B cells": "GC B cells",
    "Germinal center B cells": "GC B cells", "Age-associated B cells" : "ABCs"
}
def IGH_switched():
    
    IGH_switched = {'IGHA1': 'SW',
             'IGHA2': 'SW',
             'IGHD': 'NS',
             'IGHM': 'NS',#"#00ff7f",
             'IGHG2': 'SW',
             'IGHG3': 'SW',
             'IGHG1': 'SW',
             'IGHG4': 'SW',
             'IGHE':'SW'}
    return IGH_switched

def bcelltype_colors_dict():
   
    return colors
# Helper function for data selection
def selection_helper(df, param, value):
    selector = df[param].value_counts() > value
    idxs = selector[selector == True].index
    return(df[df[param].isin(idxs)])

# Define colors at the top level
dark_purple = "#6F4B70"
grey_green= "#ABB8B8"
sky_blue = "#78BACC"
green1 = "#77BA99"
green2 = "#68918F"
green3 = "#266967"
thistle = "#DAC1E0"
lilac = "#AF7DBB"
azul = "#0070BB"
uranian_blue = "#AFDBF5"
bice_blue = "#2072AF"
hex_colors = [dark_purple, thistle, grey_green, green1, green2, green3]
rgb_colors = [mpl.colors.to_rgb(x) for x in hex_colors]
cdict = {'red': [x[0] for x in rgb_colors], 'green': [x[1] for x in rgb_colors], 'blue': [x[2] for x in rgb_colors]}
my_cmap = mpl.colors.LinearSegmentedColormap.from_list('donor_cmap',hex_colors,N=6)
donor_cycle = [my_cmap(i/6) for i in range(6)]
donor_colors = {f'TBd{i+1}':donor_cycle[i] for i in range(len(donor_cycle))}

# Create an instance of the ColorDictionaries class
color_dicts = ColorDictionaries()
# now you can access any color dictionary by using color_dicts.get_colors(dict_name)
### IGH colors
def get_IGH_colors(simplify=False):
    if simplify:
        new_dict = {k: v for k, v in IGH_colors.items() if k[1] in ['D','M','E','s']}
        new_dict.update({'IGHA1':IGH_colors['IGHA1'],
                         'IGHA2':IGH_colors['IGHA1']}) 
        new_dict.update({'IGHG1':IGH_colors['IGHG1'],
                         'IGHG2':IGH_colors['IGHG1'],
                         'IGHG3':IGH_colors['IGHG1'],
                         'IGHG4':IGH_colors['IGHG1']})
        return new_dict
    return IGH_colors

def set_IGH_colors(new_dict):
    IGH_colors.update(new_dict)
    return IGH_colors

### donor colors
def get_donor_colors():
    return donor_colors

def set_donor_colors(new_dict):
    donor_colors.update(new_dict)
    return donor_colors

### tissue colors
def get_tissue_colors():
    return tissue_colors

def set_tissue_colors(new_dict):
    tissue_colors.update(new_dict)
    return tissue_colors

def get_celltypist_simpler():
    return celltypist_simpler


#### color palettes #####

donor_colors = {f'TBd{i+1}':donor_cycle[i] for i in range(6)}

tissue_colors = {'SP':'#801630FF', #'#765760', 
                 'PB':'#C1717B', 
                 'BM': '#E05353', #'#266967', 
                 'LN':'#E4B363',
                 'both':'k', 
                 'multiple':'k'}

IGH_colors = {'IGHA1': '#fb9a99',
             'IGHA2': '#e31a1c',
             'IGHD': '#fec44f',
             'IGHM': '#33a02c',
             'IGHG2': '#02818a',
             'IGHG4': '#6a51a3',
             'IGHG1': '#386cb0',
             'IGHG3': '#3690c0',
             'IGHE':'#67001f',
             "ambiguous":'#f7f7f7'}


bcelltype_colors_alt = {
                    'Plasmablasts':'#a6cee3',
                    'Plasma cells':'#1f78b4',
                    'Naive B cells':'#b2df8a',
                    'GC B cells':'#33a02c',
                    'ABCs':'#fb9a99',
                    'Memory B cells':'#e31a1c'}

bcelltype_colors_full = {
                    'Plasmablasts':'#a6cee3',
                    'Plasma cells':'#1f78b4',
                    'Naive B cells':'#b2df8a',
                    'Germinal center B cells':'#fdbb84',
                    'Age-associated B cells':'#fb9a99',
                    'Plasmablasts':'#a6cee3',
                    'Plasma cells':'#1f78b4',
                    'Naive B cells':'#b2df8a',
                    'Proliferative germinal center B cells':'#33a02c',
                    'Age-associated B cells':'#fb9a99',
                    'Memory B cells':'#e31a1c'}

celltypist_simpler = {
    "Proliferative germinal center B cells": "GC B cells",
    "Germinal center B cells": "GC B cells", "Age-associated B cells" : "ABCs"
}
asc_colors = {'ASC-1':'#A6CEE3',
                                         'ASC-2':'#4598C4',
                                         'ASC-3':'#1B699D',
                                         'ASC-4':'#7071A3'}

memory_b_cell_colors = {
                    'NS-Memory':'#fdb863',
                     'NS-ABCs':'#e08214',
                      'NS-Memory-CD21++':'#b35806',
                    'SW-Memory':'#b2abd2',
                    'SW-ABCs':'#8073ac',
                    'SW-Memory-CD21++':'#542788'}

#### helper functions #####

### IGH colors
def get_IGH_colors(simplify=False):
    if simplify:
        new_dict = {k: v for k, v in IGH_colors.items() if k[1] in ['D','M','E','s']}
        new_dict.update({'IGHA1':IGH_colors['IGHA1'],
                         'IGHA2':IGH_colors['IGHA1']}) 
        new_dict.update({'IGHG1':IGH_colors['IGHG1'],
                         'IGHG2':IGH_colors['IGHG1'],
                         'IGHG3':IGH_colors['IGHG1'],
                         'IGHG4':IGH_colors['IGHG1']})
        return new_dict
    return IGH_colors

def set_IGH_colors(new_dict):
    IGH_colors.update(new_dict)
    return IGH_colors

### donor colors
def get_donor_colors():
    return donor_colors

def set_donor_colors(new_dict):
    donor_colors.update(new_dict)
    return donor_colors

### tissue colors
def get_tissue_colors():
    return tissue_colors

def set_tissue_colors(new_dict):
    tissue_colors.update(new_dict)
    return tissue_colors

### cell type colors

def get_bcelltype_colors():
    return bcelltype_colors
def get_bcelltype_colors_alt():
    return bcelltype_colors_alt
def set_bcelltype__colors(new_dict):
    bcelltype_colors.update(new_dict)
    return bcelltype_colors
def get_asc_colors():
    return asc_colors
def get_mb_colors():
    return memory_b_cell_colors
def get_celltypist_simpler():
    return celltypist_simpler

def get_bcelltype_colors_full():
    return bcelltype_colors_full

def get_bcelltype_colors_alt():
    return bcelltype_colors_alt


def convert_boolean_columns(df):

    """

    Convert columns in a DataFrame that are likely to be boolean strings to boolean type.

    

    Parameters:

    - df: pd.DataFrame

        The DataFrame to process.

    

    Returns:

    - pd.DataFrame

        The DataFrame with appropriate columns converted to boolean type.

    """

    

    for col in df.columns:

        unique_values = df[col].dropna().unique()

        

        # If the unique non-null values in the column are "True" and "False", assume it's a boolean column

        if set(unique_values) == {'True', 'False'}:

            df.loc[:,col] = df[col].replace({'True': True, 'False': False}).astype('boolean')

        

        # If the column contains "True", "False", and "nan", assume it's a boolean column with NaNs

        elif set(unique_values) == {'True', 'False'} or set(unique_values) == {'True', 'False', 'nan'}:

            df.loc[:,col] = df[col].replace({'True': True, 'False': False, 'nan': np.nan}).astype('boolean')

            

    return df