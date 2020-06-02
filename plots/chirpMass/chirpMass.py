'''

\mathcal{M} =  (m_1 . m_2)^{3/5} / (m_1 + m_2)^{1/5}

'''

import numpy as np
import matplotlib.pyplot as plt

m1_length = 100
m2_length = 100

m1        =  np.ones(m1_length)
m2        =  np.ones(m2_length)
chirpMass =  np.ones([m1_length,m2_length])


for ii in range(m1_length):
    #    m1[ii] = ii+1
    m1[ii] = (ii+1) * 1e6
    
    for jj in range(m2_length):
        m2[jj] = jj+1

        #print(ii, jj, chirpMass[ii,jj])        
        chirpMass[ii,jj] = ((m1[ii] * m2[jj])**(3./5)) / ((m1[ii]+m2[jj])**(1./5))
        print(ii, m1[ii], jj, m2[jj], chirpMass[ii,jj])        



## Setting up the plot
fig, ax = plt.subplots(figsize=(8, 8), dpi=80, facecolor='w', edgecolor='k')

## Adjusting the Whitespace for the plots
left   = 0.14   # the left side of the subplots of the figure
right  = 0.98   # the right side of the subplots of the figure
bottom = 0.12   # the bottom of the subplots of the figure
top    = 0.96   # the top of the subplots of the figure
wspace = 0.16   # the amount of width reserved for blank space between subplots
hspace = 0.06   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
 
## Some NPR defaults
alpha           = 1.0
fontsize        = 16
labelsize       = fontsize
tickwidth       = 2.0
linewidth       = 2.4
tickwidth       = 2.0
ticklength      = 6.0
ticklabelsize   = labelsize
majorticklength = 12
minorticklength = 6


## `surface contour plot'
#ax.plot(m1, m2,
contours = plt.contour(m1, m2, chirpMass, 3, colors='black')
plt.clabel(contours, inline=True, fontsize=fontsize/1.4)

#plt.imshow(chirpMass, extent=[0, 100, 0, 100], origin='lower',           cmap='RdGy', alpha=0.5)
plt.imshow(chirpMass, origin='lower',           cmap='RdGy', alpha=0.5)
plt.colorbar();



## Axes limits
xmin =  0
xmax =  100
ymin =  0.0
ymax =  100
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])

## Axes labels
ax.set_xlabel(r'm_$1$ / M$_{\odot}$')
ax.set_ylabel(r'm_$2$ / M$_{\odot}$')

## Axes style
ax.tick_params(axis='both', which='major', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength,   width=tickwidth)
ax.tick_params(axis='both', which='minor', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength/2, width=tickwidth)


plt.savefig('ChirpMass_temp.png', format='png')
plt.close(fig)
        
