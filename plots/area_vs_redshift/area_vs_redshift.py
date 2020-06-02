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

## Setting the cosmology
cosmo = FlatLambdaCDM(H0=70.0, Om0=0.300)  # for round numbers!!

## So...!!! The critial bit, is that you are 'actually' passing the
## redshifts you want on the top axis to the cosmo. function, and this
## returns the DLs, i.e. the x-axis positions of ## the lables you can
## then use below.  .tolist() needed for iteration purposes.
DLs_for_axis   = np.array([0.2, 0.4, 0.6, 0.8, 1.0]) 
redshift_ticks = (cosmo.luminosity_distance(DLs_for_axis))
redshift_ticks_xx = redshift_ticks.value.tolist()

##
## LIGO 03 data 
##
path     = '/cos_pc19a_npr/data/LIGO/O3/'
filename = 'SuperEvents_full03.tbl'
table    = path+filename
LIGO_O3  = ascii.read(table)

##
## Making the plot(s)
## 
##   Distance vs.  area(50/90)
##
plt.rcParams.update({'font.size': 14})


## Some plotting defaults
ls              = 'solid'
ms              = 5
lw              = 4.0
ticklength      = 18
tickwidth       = 2.0
pointsize       = 160
pointsize_large = pointsize*1.8
fontsize        = 28
markeredgewidth = 1.5

##
##  Plotting the LIGO O3 points
##
#cmap = plt.cm.viridis
cmap = plt.cm.rainbow
#cmap = plt.cm.gist_rainbow

fig, ax = plt.subplots(figsize=(16, 8.), num=None, dpi=80, facecolor='w', edgecolor='k')
plt.style.use('dark_background')

##
##    A R E A    5 0
##
## LIGO O3 detection
ax.errorbar(         LIGO_O3['dist_Mpc'], LIGO_O3['area50'], xerr=LIGO_O3['err_dist_Mpc'],  fmt='.', color='w',zorder=-1)
hb_LIGO = ax.scatter(LIGO_O3['dist_Mpc'], LIGO_O3['area50'], s=pointsize_large, c='w')
#hb_LIGO = ax.scatter(LIGO_O3['dist_Mpc'], LIGO_O3['area50'], s=pointsize,       c=LIGO_O3['log10_FAR_Hz'],cmap=cmap)
hb_LIGO = ax.scatter(LIGO_O3['dist_Mpc'], LIGO_O3['area50'], s=pointsize,       c=LIGO_O3['BBH'],cmap=cmap)

## Colorbar and label
cbaxes = fig.add_axes([0.62, 0.28, 0.24, 0.025])
cb = fig.colorbar(hb_LIGO,  ax=ax, cax=cbaxes, orientation='horizontal') #, ticklocation = 'top')
#ax.text(5200., 2.5, 'log10(FAR/Hz)', size=fontsize/1.4, color='w')
ax.text(5500., 2.5, 'BBH', size=fontsize/1.4, color='w')

## Axes ranges; logs
#xmin = 0.00; xmax = 7300.0;  ymin = (LIGO_O3['area50'].min()/2.2);  ymax = (LIGO_O3['area50'].max()*2.2)
xmin = 0.00; xmax = 7300.0;  ymin = 1.0;  ymax = (LIGO_O3['area50'].max()*2.2)
ax.set_yscale('log')
ax.axis([xmin, xmax, ymin, ymax])

#ax.yaxis.tick_right()
ax.yaxis.set_ticks_position('both')

# Both top and bottom x-axis (must) have the same range to make sense.
zmin = xmin
zmax = xmax
ax2 = ax.twiny()
ax2.set_xticks(redshift_ticks_xx)
ax2.set_xticklabels(['{:.2}'.format(DL) for DL in DLs_for_axis], fontsize=fontsize/1.15)
ax2.set_xlim(zmin, zmax)
ax2.set_xlabel(r"redshift (H$_{0}$=70.0)", fontsize=fontsize)

ax.set_xlabel(r"Luminosity Distance / Mpc", fontsize=fontsize)
ax.set_ylabel(r"area50 / deg$^2$",          fontsize=fontsize)
plt.savefig('area50_vs_redshift_temp.png', format='png')
#plt.show()
plt.close()


##
##   a r e a   9 0
##
fig, ax = plt.subplots(figsize=(16, 8.), num=None, dpi=80, facecolor='w', edgecolor='k')
plt.style.use('dark_background')

ax.errorbar(         LIGO_O3['dist_Mpc'], LIGO_O3['area90'], xerr=LIGO_O3['err_dist_Mpc'], fmt='.', color='w',zorder=-1)
hb_LIG0 = ax.scatter(LIGO_O3['dist_Mpc'], LIGO_O3['area90'],                               s=pointsize_large, c='w')
hb_LIG0 = ax.scatter(LIGO_O3['dist_Mpc'], LIGO_O3['area90'], c=LIGO_O3['log10_FAR_Hz'],    s=pointsize,       cmap=cmap)
    

cbaxes = fig.add_axes([0.62, 0.28, 0.24, 0.025])
cb = fig.colorbar(hb_LIG0,  ax=ax, cax=cbaxes, orientation='horizontal') #, ticklocation = 'top')
#cb.set_label('redshift            ', labelpad=14)
ax.text(5200., 20., 'log10(FAR/Hz)', size=fontsize/1.4, color='w')
#ax.text(5500., 25., 'BBH', size=fontsize/1.4, color='w')

ymin = 10.
ymax = 41252. ## may as well do full sky since S190910h has A90=24264
ax.axis([xmin, xmax, ymin, ymax])
ax.set_yscale('log')
ax.yaxis.set_ticks_position('both')


# Both top and bottom x-axis (must) have the same range to make sense.
zmin = xmin
zmax = xmax
ax2 = ax.twiny()
ax2.set_xticks(redshift_ticks_xx)
ax2.set_xticklabels(['{:.2}'.format(DL) for DL in DLs_for_axis], fontsize=fontsize/1.15)
ax2.set_xlim(zmin, zmax)
ax2.set_xlabel(r"redshift (H$_{0}$=70.0)", fontsize=fontsize)

ax.set_xlabel(r"Luminosity Distance / Mpc", fontsize=fontsize)
ax.set_ylabel(r"area90 / deg$^2$",          fontsize=fontsize)
plt.savefig('area90_vs_redshift_temp.png', format='png')
plt.close()
