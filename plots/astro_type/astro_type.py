"""
=============================================
Discrete distribution as horizontal bar chart
=============================================

Stacked bar charts can be used to visualize discrete distributions.

This example visualizes the result of a survey in which people could rate
their agreement to questions on a five-element scale.

The horizontal stacking is achieved by calling `~.Axes.barh()` for each
category and passing the starting point as the cumulative sum of the
already drawn bars via the parameter ``left``.
"""

import numpy as np
import matplotlib.pyplot as plt

#category_names = ['Binary NS', 'NS-BH', 'Binary BH', 'Mass Gap', 'Terrestrial']
category_names = ['BNS', 'NS-BH', 'BBH', 'MassGap', 'Terr']

## Events by date
results_byDate = {
        '190408 ': [0.0,    0.0,	1.0,	0.0,	0.0], 
        '190412 ': [0.0,    0.0,	1.0,	0.0,	0.0], 
        '190421 ': [0.0,    0.0,	0.9674, 0.0,    0.0326],  
        '190425 ': [0.9994, 0.0,    0.0,    0.0,    0.0006], 
        '190426 ': [0.4932, 0.1293, 0.0,	0.2374, 0.1401], 
        '190503 ': [0.0,    0.0047, 0.9628, 0.0323, 0.0001],
        '190510 ': [0.9797, 0.0,    0.0,    0.0,    0.0203], 
        '190512 ': [0.0,    0.0,    0.9899, 0.0,    0.0101], 
        '190513 ': [0.0,    0.0052, 0.9427, 0.0521, 0.0],
        '190517 ': [0.0,    0.0008, 0.9826, 0.0166, 0.0],
        '190519 ': [0.0,    0.0,	0.9558, 0.0,	0.0442],
        '190521g': [0.0,    0.0,	0.9659, 0.0,	0.0341],
        '190521r': [0.0,    0.0,	0.9993, 0.0,	0.0007],
        '190602 ': [0.0,    0.0,	0.9903, 0.0,	0.0097],
        '190630 ': [0.0,    0.0052, 0.9427, 0.0521, 0.0],
        '190701 ': [0.0,    0.0,	0.9344, 0.0,	0.0656],
        '190706 ': [0.0,    0.0,	0.9898, 0.0,	0.0102],
        '190707 ': [0.0,	0.0,    1.0,    0.0,    0.0], 
        '190718 ': [0.0207, 0.0,	0.0,	0.0,	0.9793], 
        '190720 ': [0.0,    0.0,	0.9893, 0.0,	0.0107],
        '190727 ': [0.0,	0.0019,	0.9221,	0.0282,	0.0478],
        '190728 ': [0.0,	0.1435,	0.3400,	0.5165,	0.0],
        '190814 ': [0.0,   	0.9979,	0.0,	0.0021,	0.0],
        '190828 ': [0.0, 	0.0,	1.0, 	0.0, 	0.0],
        '190828 ': [0.0,	0.0,	1.0,	0.0,	0.0],
        '190901 ': [0.8607,	0.0,	0.0,	0.0,	0.1393],
        '190910d': [0.0,	0.9759,	0.0,	0.0,	0.0241],
        '190910h': [0.6122,	0.0,	0.0,	0.0,	0.3878],
        '190915 ': [0.0,	0.0,	0.9947,	0.0,	0.0053],
        '190923 ': [0.0,	0.6778,	0.0,	0.0,	0.3222],
        '190924 ': [0.0,    0.0,	0.0,	1.0,	0.0],
        '190930s': [0.0,	0.0,	0.0,	0.9508,	0.0492],
        '190930t': [0.0,	0.7426, 0.0,	0.0,	0.2574],
        '191105 ': [0.0,	0.0,	0.9531,	0.0,	0.0469],
        '191109 ': [0.0,	0.0,	1.0,	0.0,	0.0],
        '191129 ': [0.0,	0.0,	1.0,	0.0,	0.0],
        '191204':  [0.0,	0.0,	1.0,	0.0,	0.0],
        '191205 ': [0.0,	0.9321,	0.0,	0.0,	0.0679],
        '191213 ': [0.7680, 0.0,	0.0,	0.0,	0.2320],
        '191215 ': [0.0,	0.0,	0.9972,	0.0,	0.0028],
        '191216 ': [0.0,	0.0,	0.9907,	0.0093,	0.0],						
        '191222 ': [0.0,	0.0,	0.9999,	0.0,	0.0001],
        '200105 ': [0.0,	0.0273,	0.0,	0.0,	0.9727],
        '200112 ': [0.0,	0.0,	0.9997,	0.0,	0.0003],
        '200115 ': [0.0,	0.0,	0.0,	0.9997,	0.0003],
        '200128 ': [0.0,	0.0,	0.9690, 0.0,	0.0310],
        '200129 ': [0.0,	0.0,	1.00,	0.0,	0.0],
        '200208 ': [0.0,	0.0,	0.9934,	0.0,	0.0066],
        '200213 ': [0.6295, 0.0,	0.0,	0.0,	0.3705],
        '200219 ': [0.0,	0.0,	0.9640,	0.0,	0.0360],
        '200224 ': [0.0,	0.0,	0.9999, 0.0,	0.0001],
        '200225 ': [0.0,	0.0,	0.9566,	0.0,	0.0434],
        '200302 ': [0.0,	0.0,	0.8896,	0.0,	0.1104],
        '200311 ': [0.0,	0.0,	1.00,	0.0,	0.0],	
        '200316 ': [0.0,	0.0,	0.0,	0.9957,	0.0042],
}

    
## Events by type; BBH  BHS, NS-BH, MassGap, Terr.
results_byType = {
        '190408 ': [0.0,    0.0,	1.0,	0.0,	0.0], 
        '190412 ': [0.0,    0.0,	1.0,	0.0,	0.0],
        '190707 ': [0.0,	0.0,    1.0,    0.0,    0.0],
        '190828j': [0.0, 	0.0,	1.0, 	0.0, 	0.0],
        '190828l': [0.0,	0.0,	1.0,	0.0,	0.0],
        '191109 ': [0.0,	0.0,	1.0,	0.0,	0.0],
        '191129 ': [0.0,	0.0,	1.0,	0.0,	0.0],
        '191204':  [0.0,	0.0,	1.0,	0.0,	0.0],
        '200129 ': [0.0,	0.0,	1.00,	0.0,	0.0],
        '200311 ': [0.0,	0.0,	1.00,	0.0,	0.0],
        '191222 ': [0.0,	0.0,	0.9999,	0.0,	0.0001],
        '200224 ': [0.0,	0.0,	0.9999, 0.0,	0.0001],
        '200112 ': [0.0,	0.0,	0.9997,	0.0,	0.0003],
        '190521r': [0.0,    0.0,	0.9993, 0.0,	0.0007],
        '191215 ': [0.0,	0.0,	0.9972,	0.0,	0.0028],
        '190915 ': [0.0,	0.0,	0.9947,	0.0,	0.0053],
        '200208 ': [0.0,	0.0,	0.9934,	0.0,	0.0066],
        '191216 ': [0.0,	0.0,	0.9907,	0.0093,	0.0],
        '190602 ': [0.0,    0.0,	0.9903, 0.0,	0.0097],
        '190512 ': [0.0,    0.0,    0.9899, 0.0,    0.0101],
        '190706 ': [0.0,    0.0,	0.9898, 0.0,	0.0102],
        '190720 ': [0.0,    0.0,	0.9893, 0.0,	0.0107],
        '190517 ': [0.0,    0.0008, 0.9826, 0.0166, 0.0],
        '200128 ': [0.0,	0.0,	0.9690, 0.0,	0.0310],
        '190421 ': [0.0,    0.0,	0.9674, 0.0,    0.0326],
        '200219 ': [0.0,	0.0,	0.9640,	0.0,	0.0360],
        '190503 ': [0.0,    0.0047, 0.9628, 0.0323, 0.0001],
        '190519 ': [0.0,    0.0,	0.9558, 0.0,	0.0442],
        '200225 ': [0.0,	0.0,	0.9566,	0.0,	0.0434],
        '190521g': [0.0,    0.0,	0.9659, 0.0,	0.0341],
        '191105 ': [0.0,	0.0,	0.9531,	0.0,	0.0469],
        '190630 ': [0.0,    0.0052, 0.9427, 0.0521, 0.0],
        '190513 ': [0.0,    0.0052, 0.9427, 0.0521, 0.0],
        '190701 ': [0.0,    0.0,	0.9344, 0.0,	0.0656],
        '190727 ': [0.0,	0.0019,	0.9221,	0.0282,	0.0478],
        '200302 ': [0.0,	0.0,	0.8896,	0.0,	0.1104],
        '190425 ': [0.9994, 0.0,    0.0,    0.0,    0.0006],
        '190510 ': [0.9797, 0.0,    0.0,    0.0,    0.0203], 
        '190901 ': [0.8607,	0.0,	0.0,	0.0,	0.1393],
        '191213 ': [0.7680, 0.0,	0.0,	0.0,	0.2320],
        '200213 ': [0.6295, 0.0,	0.0,	0.0,	0.3705],
        '190910h': [0.6122,	0.0,	0.0,	0.0,	0.3878],
        '190426 ': [0.4932, 0.1293, 0.0,	0.2374, 0.1401],
        '190924 ': [0.0,    0.0,	0.0,	1.0,	0.0],
        '200316 ': [0.0,	0.0,	0.0,	0.9957,	0.0042],
        '200115 ': [0.0,	0.0,	0.0,	0.9997,	0.0003],
        '190930s': [0.0,	0.0,	0.0,	0.9508,	0.0492],
        '190814 ': [0.0,   	0.9979,	0.0,	0.0021,	0.0],
        '190910d': [0.0,	0.9759,	0.0,	0.0,	0.0241],
        '191205 ': [0.0,	0.9321,	0.0,	0.0,	0.0679],
        '190930t': [0.0,	0.7426, 0.0,	0.0,	0.2574],
        '190923 ': [0.0,	0.6778,	0.0,	0.0,	0.3222],
        '190728 ': [0.0,	0.1435,	0.3400,	0.5165,	0.0],
        '190718 ': [0.0207, 0.0,	0.0,	0.0,	0.9793], 
        '200105 ': [0.0,	0.0273,	0.0,	0.0,	0.9727],
}

    
#results = results_byDate
results = results_byType


def survey(results, category_names):
    """
    Parameters
    ----------
    results : dict
        A mapping from question labels to a list of answers per category.
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_names*.
    category_names : list of str
        The category labels.
    """
    plt.style.use('dark_background')

    labels = list(results.keys())
    data = np.array(list(results.values()))
    #data = np.log10(data)
    data_cum = data.cumsum(axis=1)

    ## Choosig a color-scheme 
    #category_colors = plt.get_cmap('RdYlGn_r')(np.linspace(0.15, 0.85, data.shape[1]))
    category_colors = plt.get_cmap('rainbow')(np.linspace(0.15, 0.85, data.shape[1]))
    #category_colors = plt.get_cmap('nipy_spectral')(np.linspace(0.15, 0.85, data.shape[1]))

    fig, ax = plt.subplots(figsize=(9.2, 16.9))

    ## Adjusting the Whitespace for the plots
    left   = 0.14   # the left side of the subplots of the figure
    right  = 0.94   # the right side of the subplots of the figure
    bottom = 0.04   # the bottom of the subplots of the figure
    top    = 0.96   # the top of the subplots of the figure
    wspace = 0.26   # the amount of width reserved for blank space between subplots
    hspace = 0.06   # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

    ## Some NPR defaults
    lw              = 1.0
    ls              = 'solid'
    ms              = 1.
    ms_large        = ms*8.
    alpha           = 1.0
    fontsize        = 16
    labelsize       = fontsize
    tickwidth       = 2.0
    linewidth       = 2.4
    tickwidth       = 2.0
    ticklength      = 6.0
    ticklabelsize   = labelsize
    majorticklength = 12
    minorticklength = 6


    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        ax.barh(labels, widths, left=starts, height=0.5,
                label=colname, color=color)
        xcenters = starts + widths / 2

        r, g, b, _ = color
        text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
        #for y, (x, c) in enumerate(zip(xcenters, widths)):
        #    ax.text(x, y, str(int(c)), ha='center', va='center',color=text_color)
#    ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1), loc='lower left', fontsize='x-large')
    ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1), loc='lower left', fontsize=fontsize)

    plt.rcParams['text.usetex'] = True
    plt.text(0.75, 56.0, r'{\bf @npr247}',  fontsize=fontsize*1.8)

    return fig, ax


survey(results, category_names)





#plt.show()
plt.savefig('astro_type_temp.png',format='png')
#plt.savefig('astro_type_temp.pdf',format='pdf')

plt.close()






