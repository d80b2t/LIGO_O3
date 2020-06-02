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
#cosmo = FlatLambdaCDM(H0=67.7, Om0=0.307)  # Banados thesis
cosmo = FlatLambdaCDM(H0=69.6, Om0=0.286)  # currnet Ned Wright CosmoCalc default

## Hmmm, can't get this to work right now...
#DLs_for_axis = np.array([0.01, 500., 1000., 2000., 3000., 5000., 6000.]) * u.Mpc
#DLs_for_axis   = np.array([2000., 4000., 6000.]) * u.Mpc
#redshift_ticks = [z_at_value(cosmo.luminosity_distance, DLs) for DLs in DLs_for_axis]

DLs_for_axis   = np.array([0.2, 0.4, 0.6, 0.8, 1.0]) 
redshift_ticks = (cosmo.luminosity_distance(DLs_for_axis))
redshift_ticks_xx = redshift_ticks.value.tolist()

## from my ol' Lz.py script
ages = np.array([13, 10, 8, 6, 5, 4, 3, 2, 1.5, 1.2, 1, 0.8, 0.70]) * u.Gyr
ageticks = [z_at_value(cosmo.age, age) for age in ages]


##
## LIGO 03 data 
##
path     = '/cos_pc19a_npr/data/LIGO/O3/'
filename = 'SuperEvents_bayestar.tbl'
table    = path+filename
LIGO_O3  = ascii.read(table)


##
## Making the plot(s)
## 
##  Distance vs.  area(50/90)
##
plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(figsize=(10, 8.), num=None, dpi=80, facecolor='w', edgecolor='k')

plt.style.use('dark_background')

## Some plotting defaults
ls              = 'solid'
ms              = 5
lw              = 4.0
ticklength      = 18
tickwidth       = 2.0
pointsize       = 100
pointsize_large = pointsize*2.2
fontsize        = 28
markeredgewidth = 1.5

##
##  Plotting the LIGO O3 points
##
#cmap = plt.cm.rainbow
cmap = plt.cm.gist_rainbow

ax.errorbar(LIGO_O3['distance_Mpc'], LIGO_O3['area50'], xerr=LIGO_O3['err_dist_Mpc'],  
                         #markeredgecolor ='w', markeredgewidth = markeredgewidth,
                        fmt='.', color='w',zorder=-1)

## Good W3 detections
hb_VHzQ = ax.scatter(LIGO_O3['distance_Mpc'], LIGO_O3['area50'],
                     s=pointsize_large, c='w')


hb_VHzQ = ax.scatter(LIGO_O3['distance_Mpc'], LIGO_O3['area50'], 
                     c=LIGO_O3['log10_FAR_Hz'],
                     s=pointsize, cmap=cmap)



cbaxes = fig.add_axes([0.52, 0.76, 0.34, 0.025])
cb = fig.colorbar(hb_VHzQ,  ax=ax, cax=cbaxes, orientation='horizontal') #, ticklocation = 'top')
#cb.set_label('redshift            ', labelpad=14)
ax.text(3600.00, 2000.0, 'log10(FAR/Hz)', size=fontsize/1.4, color='w')

xmin =  0.00; xmax =  7300.0; ymin = 5.1;  ymax =  4100.
ax.set_yscale('log')
ax.axis([xmin, xmax, ymin, ymax])

#ax.yaxis.tick_right()
ax.yaxis.set_ticks_position('both')


#ax.tick_params('both', direction='in', which='major', bottom=True, top=True, left=True, right=True,               length=majorticklength, width=tickwidth)

zmin = xmin
zmax = xmax
ax2 = ax.twiny()
ax2.set_xticks(redshift_ticks_xx)
ax2.set_xticklabels(['{:g}'.format(DL) for DL in DLs_for_axis], fontsize=fontsize/1.2)
ax2.set_xlim(zmin, zmax)

#ax4 = ax.twiny()
#ax4.set_xticks(ageticks)
#ax4.set_xticklabels(['{:g}'.format(age) for age in ages.value], fontsize=fontsize/1.2)


ax.set_xlabel(r"Luminosity Distance / Mpc", fontsize=fontsize)
ax.set_ylabel(r"area50 / deg$^2$", fontsize=fontsize)

plt.savefig('area50_vs_redshift_temp.png', format='png')
#plt.show()
plt.close()


#DLs_for_axis = np.array([2000., 4000., 6000.]) * u.Mpc
#redshift_ticks = [z_at_value(cosmo.luminosity_distance, DLs) for DLs in DLs_for_axis]



# ages = np.array([13, 10, 8, 6, 5, 4, 3, 2, 1.5, 1.2, 1, 0.8, 0.70]) * u.Gyr
# ageticks = [z_at_value(cosmo.age, age) for age in ages]

#ax4 = main_ax.twiny()
#ax4.set_xticks(ageticks)
#ax4.set_xticklabels(['{:g}'.format(age) for age in ages.value], fontsize=fontsize/1.2)
