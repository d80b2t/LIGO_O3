import sys
import math
import json
import numpy
import healpy
from colorsys import hsv_to_rgb
import matplotlib
#matplotlib.use('Agg') # switch to agg to suppress plotting windows 
import matplotlib.pyplot as plt


radian = 180/numpy.pi

filename = 'bayestar_S190521r.fits.gz'

healpix_data, distmu, distsigma, distnorm = healpy.read_map(filename, field=[0, 1, 2, 3])

dim = 3

print("dimension = %d" % dim)

if dim == 3:
    ii = (~numpy.isinf(distmu)) & (distmu > 0.0)
    print("good percent %.2f" % (numpy.sum(ii)*100.0/len(distmu)))
    dm = distmu[ii]
    prob = healpix_data[ii]
    x = numpy.multiply(prob,dm)
    print("Average distance %.2f " % (numpy.sum(x)/numpy.sum(prob)))

