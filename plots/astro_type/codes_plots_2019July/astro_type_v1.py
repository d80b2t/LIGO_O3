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


category_names = ['BNS', 'NSBH', 'BBH', 'MassGap', 'Terrestrial']
'''
results = {
        '190408an': [0.0,    0.0,	 1.0,	 0.0,	 0.0], 
        '190412m ':	[0.0,    0.0,	 1.0,	 0.0,	 0.0], 
        '190421ar':	[0.0,    0.0,	 0.9674, 0.0,    0.0326],  
        '190425z ': [0.9994, 0.0,    0.0,    0.0,    0.0006], 
        '190426c ':	[0.4932, 0.1293, 0.0,	 0.2374, 0.1401], 
        '190503bf':	[0.0,    0.0047, 0.9628, 0.0323, 0.0001],
        '190510g ':	[0.9797, 0.0,    0.0,    0.0,    0.0203], 
        '190512at': [0.0,    0.0,    0.9899, 0.0,    0.0101], 
        '190513bm':	[0.0,    0.0052, 0.9427, 0.0521, 0.0],
        '190517h ':	[0.0,    0.0008, 0.9826, 0.0166, 0.0],
        '190519bj':	[0.0,    0.0,	 0.9558, 0.0,	 0.0442],
        '190521g ':	[0.0,    0.0,	 0.9659, 0.0,	 0.0341],
        '190521r ':	[0.0,    0.0,	 0.9993, 0.0,	 0.0007],
        '190602aq':	[0.0,    0.0,	 0.9903, 0.0,	 0.0097],
        '190630ag':	[0.0,    0.0052, 0.9427, 0.0521, 0.0],
        '190701ah':	[0.0,    0.0,	 0.9344, 0.0,	 0.0656],
        '190706ai':	[0.0,    0.0,	 0.9898, 0.0,	 0.0102],
        '190707q ':	[0.0,	 0.0,    1.0,    0.0,    0.0], 
        '190718y ': [0.0207, 0.0,	 0.0,	 0.0,	 0.9793], 
        '190720a ':	[0.0,    0.0,	 0.9893, 0.0,	 0.0107]
}
'''
results = {
        '190408': [0.0,    0.0,	 1.0,	 0.0,	 0.0], 
        '190412':	[0.0,    0.0,	 1.0,	 0.0,	 0.0], 
        '190421':	[0.0,    0.0,	 0.9674, 0.0,    0.0326],  
        '190425': [0.9994, 0.0,    0.0,    0.0,    0.0006], 
        '190426':	[0.4932, 0.1293, 0.0,	 0.2374, 0.1401], 
        '190503':	[0.0,    0.0047, 0.9628, 0.0323, 0.0001],
        '190510':	[0.9797, 0.0,    0.0,    0.0,    0.0203], 
        '190512': [0.0,    0.0,    0.9899, 0.0,    0.0101], 
        '190513':	[0.0,    0.0052, 0.9427, 0.0521, 0.0],
        '190517':	[0.0,    0.0008, 0.9826, 0.0166, 0.0],
        '190519':	[0.0,    0.0,	 0.9558, 0.0,	 0.0442],
        '190521':	[0.0,    0.0,	 0.9659, 0.0,	 0.0341],
        '190521':	[0.0,    0.0,	 0.9993, 0.0,	 0.0007],
        '190602':	[0.0,    0.0,	 0.9903, 0.0,	 0.0097],
        '190630':	[0.0,    0.0052, 0.9427, 0.0521, 0.0],
        '190701':	[0.0,    0.0,	 0.9344, 0.0,	 0.0656],
        '190706':	[0.0,    0.0,	 0.9898, 0.0,	 0.0102],
        '190707':	[0.0,	 0.0,    1.0,    0.0,    0.0], 
        '190718': [0.0207, 0.0,	 0.0,	 0.0,	 0.9793], 
        '190720':	[0.0,    0.0,	 0.9893, 0.0,	 0.0107]
}


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
    #category_colors = plt.get_cmap('RdYlGn_r')(np.linspace(0.15, 0.85, data.shape[1]))
    #category_colors = plt.get_cmap('rainbow')(np.linspace(0.15, 0.85, data.shape[1]))
    category_colors = plt.get_cmap('plasma_r')(np.linspace(0.15, 0.85, data.shape[1]))

    fig, ax = plt.subplots(figsize=(9.2, 6.5))
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
    ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1), loc='lower left', fontsize='x-large')

    return fig, ax


survey(results, category_names)

#plt.show()
plt.savefig('astro_type_temp.png',format='png')
plt.savefig('astro_type_temp.pdf',format='pdf')

plt.close()






