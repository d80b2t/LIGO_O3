'''
Trying to combine (either the BAYESTAR or LALInference) LIGO/Virgo skymaps
Using code in:
  Singer_2016_ApJS_226_10.pdf
  
The first column, PROB, is simply the probability that the source is
contained within the pixel i that is centered on the direction n_i,
the same as in the 2D localization format.
The second and third columns, DISTMU and DISTSTD, are the ansatz location and scale parameters, respectively.
The fourth column, DISTNORM, is the ansatz normalization coefficient, included for convenience.

PROB      pixel^-1     Probability that the source is contained in pixel i, centered on the direction n_i
DISTMU    Mpc          Ansatz location parameter of conditional distance distribution in direction n_i, or \infty if invalid
DISTSIGMA Mpc          Ansatz scale parameter of conditional distance distribution in direction n_i, or 1 if invalid
DISTNORM  Mpc^-2       Ansatz normalization coefficient, or 0 if invalid


First go is goinng to be for   GW190828j and GW190828l

'''


import healpy as hp
import numpy as np

from matplotlib import pyplot as plt
from scipy.stats import norm

from ligo.skymap.postprocess import find_greedy_credible_levels


datapath = '/cos_pc19a_npr/data/LIGO/O3/'

S190828l_file = datapath+'S190828l/bayestar.fits.gz'
S190828j_file = datapath+'S190828j/bayestar.fits.gz'


S190828l_prob, S190828l_distmu, S190828l_distsigma, S190828l_distnorm = hp.read_map(S190828l_file, field = [0, 1, 2, 3])
print('S190828l/bayestar.fits.gz read-in...')

S190828j_prob, S190828j_distmu, S190828j_distsigma, S190828j_distnorm = hp.read_map(S190828j_file, field = [0, 1, 2, 3])
print('S190828j/bayestar.fits.gz read-in...')
print()
print()

#hpx, header = hp.read_map(inputdata, h=True)


S190828l_npix = len(S190828l_prob)
S190828j_npix = len(S190828j_prob)

S190828l_nside = hp.npix2nside(S190828l_npix)
S190828j_nside = hp.npix2nside(S190828j_npix)

S190828l_pixarea_deg2 = hp.nside2pixarea(S190828l_nside, degrees = True)
S190828j_pixarea_deg2 = hp.nside2pixarea(S190828j_nside, degrees = True)



S190828l_credible_levels = find_greedy_credible_levels(S190828l_prob)
S190828l_npix50    = np.count_nonzero(S190828l_credible_levels <= 0.50)
S190828l_area50    = np.sum(S190828l_credible_levels           <= 0.50) * hp.nside2pixarea(S190828l_nside, degrees=True)
S190828l_boolArr50 =       (S190828l_credible_levels           <= 0.50)
S190828l_result50  = np.where(S190828l_boolArr50)
S190828l_npix90    = np.count_nonzero(S190828l_credible_levels <= 0.90)
S190828l_area90    = np.sum(S190828l_credible_levels           <= 0.90) * hp.nside2pixarea(S190828l_nside, degrees=True)
S190828l_boolArr90 =       (S190828l_credible_levels           <= 0.90)
S190828l_result90  = np.where(S190828l_boolArr90)
print('S190828l: npix50, area50  ', S190828l_npix50, S190828l_area50,   '  npix90, area90  ', S190828l_npix90, S190828l_area90) 


S190828j_credible_levels = find_greedy_credible_levels(S190828j_prob)
S190828j_npix50    = np.count_nonzero(S190828j_credible_levels <= 0.50)
S190828j_area50    = np.sum(S190828j_credible_levels           <= 0.50) * hp.nside2pixarea(S190828j_nside, degrees=True)
S190828j_boolArr50 =       (S190828j_credible_levels           <= 0.50)
S190828j_result50  = np.where(S190828j_boolArr50)
S190828j_npix90    = np.count_nonzero(S190828j_credible_levels <= 0.90)
S190828j_area90    = np.sum(S190828j_credible_levels           <= 0.90) * hp.nside2pixarea(S190828j_nside, degrees=True)
S190828j_boolArr90 =       (S190828j_credible_levels           <= 0.90)
S190828j_result90  = np.where(S190828j_boolArr90)
print('S190828j: npix50, area50  ', S190828j_npix50, S190828j_area50,   '  npix90, area90  ', S190828j_npix90, S190828j_area90) 
print()



## The figure...
#fig, ax = plt.subplots(figsize=(5, 3), dpi=80, facecolor='w', edgecolor='k')

## Plotting a Mollweide-projection all-sky map:
hp.mollview(S190828l_prob)
plt.savefig("S190828l_prob_temp.png")

hp.mollview(np.log10(S190828j_prob))
plt.savefig("S190828l_prob_log10_temp.png")


hp.mollview((S190828j_prob))
plt.savefig("S190828j_prob_temp.png")

hp.mollview(np.log10(S190828j_prob))
plt.savefig("S190828j_prob_log10_temp.png")


hp.mollview(S190828l_prob*S190828j_prob)
plt.savefig("two_probs_product_temp.png")
hp.mollview((np.log10(S190828l_prob))+(np.log10(S190828j_prob)))
plt.savefig("two_probs_log10_sum_temp.png")


##
product_prob = S190828l_prob+S190828j_prob

product_npix = len(product_prob)
product_nside = hp.npix2nside(product_npix)

product_credible_levels = find_greedy_credible_levels(product_prob)
product_npix50 = np.count_nonzero(product_credible_levels <= 0.50)
product_area50 = np.sum(product_credible_levels           <= 0.50) * hp.nside2pixarea(product_nside, degrees=True)
product_npix90 = np.count_nonzero(product_credible_levels <= 0.90)
product_area90 = np.sum(product_credible_levels           <= 0.90) * hp.nside2pixarea(product_nside, degrees=True)

print('product: npix50, area50  ', product_npix50, product_area50,   '  npix90, area90  ', product_npix90, product_area90) 
print()


'''
S190828l >
LOGBSN  =    20.43991117673434 / Log Bayes factor: signal vs. noise             

S190828j > fitsheader bayestar.fits.gz
LOGBSN  =    66.38571899172567 / Log Bayes factor: signal vs. noise             
'''


#hp.projplot(theta, phi, 'bo')  # plot 'o' in blue at coord (theta, phi)
#fig.savefig("test.png")
#plt.show()
