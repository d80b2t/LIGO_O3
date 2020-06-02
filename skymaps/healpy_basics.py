'''
http://balbuceosastropy.blogspot.com/2013/09/representations-of-sphere-and-cosmic.html
'''

from __future__ import division

import numpy as np
import healpy as hp
import astroML

#This removes some nasty deprecation warnings that do not interfere with the execution
import warnings
warnings.filterwarnings('ignore')

'''
# This IPython magic generates a table with version information
#https://github.com/jrjohansson/version_information
%load_ext version_information
%version_information numpy, healpy, astroML, astroPy
'''

print()
for NSIDE in 2.**np.arange(12):
    print('The number of pixels for NSIDE = %2d is %5d' %(NSIDE, hp.nside2npix(NSIDE)))
print()



