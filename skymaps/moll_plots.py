'''
https://stackoverflow.com/questions/37844221/healpy-plotting-how-do-i-make-a-figure-with-subplots-using-the-healpy-mollview

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(ncols=2)

plt.axes(ax1)
hp.mollview(np.random.random(hp.nside2npix(32)), hold=True)

plt.axes(ax2)
hp.mollview(np.arange(hp.nside2npix(32)), hold=True)

'''


'''
https://stackoverflow.com/questions/57560428/aitoff-projections-using-healpy-and-projplot/57561507#57561507
'''

import healpy as hp
import numpy as np

nside = 64
npix = hp.nside2npix(nside)
arr = np.random.randn(npix)

# Draw a circle
r = np.full(100, 20.)
phi = np.linspace(0., 2*np.pi, 100)
x = np.cos(phi)*r
y = np.sin(phi)*r

# Plot the map and the circle
hp.mollview(arr)
hp.projplot(x, y, c='r', lonlat=True)

