'''
https://dcc.ligo.org/public/0119/P1500071/007/index.html

'''


import healpy as hp
import numpy as np

from matplotlib import pyplot as plt
from scipy.stats import norm
from astroquery.ned import Ned
from astropy.utils.data import download_file
from astropy.cosmology import default_cosmology

cosmo = default_cosmology.get()

url = 'https://dcc.ligo.org/P1500071/public/18951_bayestar.fits.gz'
filename = download_file(url, cache=True)

prob, distmu, distsigma, distnorm = hp.read_map(filename, field=[0, 1, 2, 3])

