'''
Make a nice plot of LIGO O3 (and maybe even O1/O2 events) 
'''
import math
import numpy as np

from astropy.io        import fits, ascii
from astropy.table     import Table
from astropy.cosmology import z_at_value

import fitsio
from fitsio import FITS,FITSHDR

import matplotlib
import matplotlib.pyplot as plt
from matplotlib        import colors as mcolors
from matplotlib        import gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator


## Setting up the cosmology
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

# In this case we just need to define the matter density 
# and hubble parameter at z=0.
# Note the default units for the hubble parameter H0 are km/s/Mpc. 
# You can also pass an astropy `Quantity` with the units specified. 
cosmo = FlatLambdaCDM(H0=67.7, Om0=0.307)  #Banados thesis

## Hmmm, can't get this to work right now...
DLs_for_axis = np.array([0.01, 500., 1000., 2000., 3000., 5000., 6000.]) * u.Mpc
redshift_ticks = [z_at_value(cosmo.luminosity_distance, DLs) for DLs in DLs_for_axis]

##
## LIGO 03 data 
##
path     = '/cos_pc19a_npr/data/LIGO/O3/'
filename = 'SuperEvents_bayestar.tbl'
table    = path+filename
LIGO_O3  = ascii.read(table)


##
##  Making the plot(s)
## 
##  distance vs.  area(50/90)
##
plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(figsize=(10, 10), num=None, dpi=80, facecolor='w', edgecolor='k')

plt.style.use('dark_background')

## Some plotting defaults
ms              = 5
lw              = 4.0
ls              = 'solid'
ticklength      = 18
tickwidth       = 2.0
pointsize       = 100
pointsize_large = pointsize*2.2
fontsize        = 28

##
##  Plotting the LIGO O3 points
##
cmap = plt.cm.seismic

#ax.errorbar(LIGO_O3['distance_Mpc'], LIGO_O3['area50'], xerr=LIGO_O3['err_dist_Mpc'],            fmt='none', linestyle=ls, ms=ms, linewidth=lw, color='w')

#ax.scatter(LIGO_O3['area50'],                  (VHzQ['unW1mag'] - VHzQ['unW2mag']),           s=pointsize_large/2.5, c='k')
#ax.scatter((VHzQ['unW2mag'] - VHzQ['w3mpro']), (VHzQ['unW1mag'] - VHzQ['unW2mag']),           c=VHzQ['redshift'],           s=pointsize/2.5, cmap=cmap)


## Good W3 detections
hb_VHzQ = ax.scatter(LIGO_O3['distance_Mpc'], LIGO_O3['area50'],
                     s=pointsize_large, c='w')


hb_VHzQ = ax.scatter(LIGO_O3['distance_Mpc'], LIGO_O3['area50'], 
                     c=LIGO_O3['log10_FAR_Hz'],
                     s=pointsize, cmap=cmap)


cbaxes = fig.add_axes([0.18, 0.32, 0.34, 0.025])
cb = fig.colorbar(hb_VHzQ,  ax=ax, cax=cbaxes, orientation='horizontal', ticklocation = 'top')
#cb.set_label('redshift            ', labelpad=14)
ax.text(0.00, -1.0, 'redshift', size=fontsize/1.2)

xmin =  0.00; xmax =  6000.0; ymin = 10.0;  ymax =  2900.
ax.set_yscale('log')
ax.axis([xmin, xmax, ymin, ymax])

#ax.tick_params(axis='both', which='major', labelsize=lw, top='on', right='on', direction='in', length=ticklength,   width=tickwidth)
#ax.tick_params(axis='both', which='minor', labelsize=lw, top='on', right='on', direction='in', length=ticklength/2, width=tickwidth)

#majorLocator_x = MultipleLocator(1000)
#minorLocator_x = MultipleLocator(200)
#majorLocator_y = MultipleLocator(1000)
#minorLocator_y = MultipleLocator(200)

#ax.xaxis.set_major_locator(majorLocator_x)
#ax.xaxis.set_minor_locator(minorLocator_x)
#ax.yaxis.set_major_locator(majorLocator_y)
#ax.yaxis.set_minor_locator(minorLocator_y)


ax.set_xlabel(r"luminosity distance", fontsize=fontsize)
ax.set_ylabel(r"area50", fontsize=fontsize)

plt.savefig('area_vs_redshift_temp.png', format='png')
#plt.show()
plt.close()
