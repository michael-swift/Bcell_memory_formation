import matplotlib as mpl

#### hex codes #####
dark_purple = "#6F4B70"
grey_green= "#abb8b8"
sky_blue = "#78BACC"
green1 = "#77BA99"
green2 = "#68918f"
green3 = "#266967"
thistle = "#DAC1E0"

hex_colors = [dark_purple, thistle, grey_green, green1, green2, green3] 
rgb_colors = [mpl.colors.to_rgb(x) for x in hex_colors]
cdict = {'red': [x[0] for x in rgb_colors],
         'green': [x[1] for x in rgb_colors],
         'blue': [x[2] for x in rgb_colors]
         }
donor_cmap = mpl.colors.LinearSegmentedColormap.from_list('donor_cmap',hex_colors,N=6)
donor_cycle = [donor_cmap(i/6) for i in range(6)]

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



bcelltype_colors = {
                    'Plasmablasts':'#a6cee3',
                    'Plasma cells':'#1f78b4',
                    'Naive B cells':'#b2df8a',
                    'Proliferative germinal center B cells':'#33a02c',
                    'Age-associated B cells':'#fb9a99',
                    'Memory B cells':'#e31a1c',
                    }
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
                    'ASCs': '#1f78b4'}
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
    return donor_colors

def set_tissue_colors(new_dict):
    tissue_colors.update(new_dict)
    return tissue_colors

### tissue colors

def get_bcelltype_colors():
    return bcelltype_colors

def set_bcelltype_colors(new_dict):
    bcelltype_colors.update(new_dict)
    return bcelltype_colors