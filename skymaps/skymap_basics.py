'''
Huge h/t to
https://emfollow.docs.ligo.org/userguide/tutorial/skymaps.html
'''

import healpy as hp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from ligo.skymap.postprocess import find_greedy_credible_levels


path     = '/cos_pc19a_npr/data/LIGO/O3/S190814bv/'
filename = 'LALInference.v1.fits.gz'
inputdata = path+filename

print()
print('Reading in ', inputdata)
print()


hpx, header = hp.read_map(inputdata, h=True)
## Each entry in the array represents the probability contained within
## a quadrilateral pixel whose position on the sky is uniquely
## specified by the index in the array and the array’s length.

## The image data is a 1D array of values:
hpx

## Most Probable Sky Location
## Let’s find the highest probability pixel.
ipix_max = np.argmax(hpx)


## Because HEALPix pixels are equal area, we can find the number of
## pixels per square degree just from the length of the HEALPix array:
npix = len(hpx)
sky_area_deg   = 4 * 180**2 / np.pi
sky_area_arcmin = 4 * 180**2 / np.pi *(60*60.)
sky_area_deg / npix
print()
print('The sky area per pixel is ', sky_area_deg    / npix, ' deg^2' )
print('The sky area per pixel is ', sky_area_arcmin / npix, ' arcmin^2', )
print()

## https://healpix.jpl.nasa.gov
## 3,145,728 pixels (~7 arcmin resolution) and the second (right)
## is a model of the Cosmic Microwave Background (CMB) radiation
## temperature anisotropy, composed of 12,582,912 pixels (~3.4 arcmin
## resolution).

## The function hp.pix2ang converts from pixel index to spherical polar coordinates;
## The function hp.ang2pix does the reverse.

## Both hp.pix2ang and hp.ang2pix take, as their first argument,
## nside, the lateral resolution of the HEALPix map.  You can find nside
## from the length of the image array by calling hp.npix2nside:
nside = hp.npix2nside(npix)
print('nside', nside)

ra = 194.95
dec = 27.98
theta = 0.5 * np.pi - np.deg2rad(dec)
phi = np.deg2rad(ra)
ipix = hp.ang2pix(nside, theta, phi)
ipix


## Test if a Sky Location is in the 90% Credible Region
i = np.flipud(np.argsort(hpx))
sorted_credible_levels = np.cumsum(hpx[i])
credible_levels = np.empty_like(sorted_credible_levels)
credible_levels[i] = sorted_credible_levels
credible_levels

## N.B. !!
## Observe that the values in the resulting credible level map vary
## inversely with probability density: the most probable pixel is
## assigned to the credible level 0.0, and the least likely pixel is
## assigned the credible level 1.0.

credible_levels = find_greedy_credible_levels(hpx)
credible_levels

## To check if the pixel that we identified in the previous section is
## within the 90% credible level, simply test if the value of the
## credible level map is less than or equal to 0.9 at that pixel:
credible_levels[ipix]
credible_levels[ipix] <=0.9

##Find the Area of the 90% Credible Region
## Since we just found the credible level map, it’s easy to compute
## the 90% credible area by counting the number of pixels inside the
## 90% credible region and multiplying by the area per pixel.

## In the Python expression below, note that (credible_levels <= 0.9)
## evaluates to a binary array; when it is summed over, true values
## are treated as 1 and false values are treated as 0.

npix10 = np.count_nonzero(credible_levels<=.1)
npix50 = np.count_nonzero(credible_levels<=.5)
npix90 = np.count_nonzero(credible_levels<=.9)

area10 = np.sum(credible_levels <= 0.10) * hp.nside2pixarea(nside, degrees=True)
area50 = np.sum(credible_levels <= 0.50) * hp.nside2pixarea(nside, degrees=True)
area90 = np.sum(credible_levels <= 0.90) * hp.nside2pixarea(nside, degrees=True)

print()
print('npix10, npix50, npix90', npix10, npix50, npix90)
print('area10, area50, area90', area10, area50, area90)
print()


ipix = hpx[(credible_levels<=.1)]


## What are the R.A.'s and Decl's of the pixels with the e.g. area10, area50, area90 values??
## Pick out the True pixesl
boolArr10 = (credible_levels <=.1)
result10 = np.where(boolArr10)
boolArr50 = (credible_levels <=.5)
result50 = np.where(boolArr50)
boolArr90 = (credible_levels <=.9)
result90 = np.where(boolArr90)


theta, phi = hp.pix2ang(nside, result50[0])
ra = np.rad2deg(phi)
dec = np.rad2deg(0.5 * np.pi - theta)
print(ra, dec)




## https://healpy.readthedocs.io/en/latest/healpy_visu.html#tracing-lines-or-points
fig, ax = plt.subplots(figsize=(5, 3), dpi=80, facecolor='w', edgecolor='k')
#fig, ax = plt.subplots()
## Plotting a Mollweide-projection all-sky map:
hp.mollview(hpx)
hp.projplot(theta, phi, 'bo')  # plot 'o' in blue at coord (theta, phi)
fig.savefig("test.png")
plt.show()


#fig, (ax1, ax2) = plt.subplots(ncols=2)
#plt.axes(ax1)
#hp.mollview(np.random.random(hp.nside2npix(32)), hold=True)
#plt.axes(ax2)
#hp.mollview(np.arange(hp.nside2npix(32)), hold=True)
