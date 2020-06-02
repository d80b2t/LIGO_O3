'''
Straight from::
    https://lscsoft.docs.ligo.org/ligo.skymap/ligo/skymap/plot/pp.html
'''



import ligo.skymap.plot
from matplotlib import pyplot as plt
import numpy as np

n = 100
p_values_1 = np.random.uniform(size=n) # One experiment
p_values_2 = np.random.uniform(size=n) # Another experiment
p_values_3 = np.random.uniform(size=n) # Yet another experiment

fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111, projection='pp_plot')
ax.add_confidence_band(n, alpha=0.95) # Add 95% confidence band
ax.add_diagonal() # Add diagonal line
ax.add_lightning(n, 20) # Add some random realizations of n samples
ax.add_series(p_values_1, p_values_2, p_values_3) # Add our data

plt.savefig('pp_temp.png') 
