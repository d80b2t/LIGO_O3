'''
arXiv:1605.04242v3
'''

import math
import healpy as hp
import numpy as np
from matplotlib import pyplot as plt 
from scipy.stats import norm

from astropy.utils.data import download_file


## url = ('https://dcc.ligo.org/P1500071/public/18951_bayestar.fits.gz')
## filename = download_file(url, cache=True)


hpx = hp.read_map('bayestar.fits.gz,0')


prob = hp.read_map(filename)

#map,header = hp.read_map(filename, h=True)

##prob, distmu, distsigma, distnorm = hp.read_map(filename, field=[0, 1, 2, 3])
prob, distmu, distsigma, distnorm = hp.read_map('bayestar.fits.gz,0', field=[0, 1, 2, 3])


npix = len(prob)
npix

nside = hp.npix2nside(npix)
nside

ra, dec = 137.8, -39.9

theta = 0.5 * np.pi - np.deg2rad(dec)
phi = np.deg2rad(ra)


ipix = hp.ang2pix(nside, theta, phi)


r = np.linspace(0, 6000)

dp_dr = r**2 * distnorm[ipix] * norm(distmu[ipix], distsigma[ipix]).pdf(r)
plt.plot(r, dp_dr)


