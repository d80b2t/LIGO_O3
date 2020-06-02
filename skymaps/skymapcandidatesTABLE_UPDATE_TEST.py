import mysql.connector
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np
import requests
import json
from collections import OrderedDict
import os
from astropy.time import Time
from astropy.table import Table
from astropy.io import ascii
from datetime import date
os.environ["GIT_PYTHON_REFRESH"] = 'quiet'
import mosfit
import healpy as hp
import settings

msl = mysql.connector.connect(            user    =settings.DB_USER,             password=settings.DB_PASS,             host    =settings.DB_HOST,             database='ztf')


"""
Determines what the current jd time is and user inputs cut off time for data collection
"""
# dates - 08042019-15042019 most optimal
startdate = 2458581.500000 #float(input("Start Date (JD): "))
enddate = 2458588.500000 #float(input("End Date (JD): "))
# quality of candidates
candquality = 'c.rb >= 0.5 and c.elong <= 1.2 and abs(c.magdiff) <= 0.3 AND c.isdiffpos = \'t\''

"""Grab's skymap from LASAIR"""
get_ipython().system('wget https://lasair.roe.ac.uk/lasair/static/ztf/skymap/GW170817.fits.gz')
print("Done.")


"""
Identifies required data from imported skymap and finds useful pixels
"""
# resolution of skymap
NSIDE = 1024 #12 * 1024^2
# number of pixels in image
NPIX =hp.nside2npix(NSIDE)
# ordered list by pixel id
m = np.arange(NPIX)

filename = 'GW170817.fits.gz'

# assigns 4 variables with appropriate data from skymap
hpdata, distmu, distsigma, distnorm = hp.read_map(filename, field=[0, 1, 2, 3])

# determines theta and phi values for each pixel
theta, phi = np.degrees(hp.pix2ang(NSIDE,ipix = m))

# finds pixels with probability >1% of maximum probability
reg_pix_theta = theta[hpdata>max(hpdata)/100]
reg_pix_phi = phi[hpdata>max(hpdata)/100]
reg_hpdata = hpdata[hpdata>max(hpdata)/100]
reg_distmu = distmu[hpdata>max(hpdata)/100]
reg_distsigma = distsigma[hpdata>max(hpdata)/100]

print("Done.")


"""
Find's skymap data range and ensures format is suitable for SQL queries
"""
# Co-ordinates for set patch of sky
declend = float(-max(reg_pix_theta)+90)
declst = float(-min(reg_pix_theta)+90)
raend = float(max(reg_pix_phi))
rast = float(min(reg_pix_phi))
# Determines start and end points for ra coordinates
ra1 = rast
ra2 = raend
if ra1 < ra2:
    pass
else:
    ra1 = raend
    ra2 = rast
#determines start and end points for dec coordinates
decl1 = declst
decl2 = declend
if decl1 < decl2:
    pass
else:
    decl1 = declend
    decl2 = declst
print("ra corrected is from ", ra1, " to ", ra2)
print("declination corrected is from ", decl1, " to ", decl2)

mycursor = msl.cursor(buffered = True, dictionary = True)

"""Finds all transients within ra, dec, boundries where earliest object detection can be before and after cut-off. """
#bad candidate query
badquerycand = 'SELECT c.candid, c.jd, c.ra, c.decl, c.magpsf, c.fid FROM candidates c NATURAL JOIN objects o WHERE '
badquerycand += 'c.ra BETWEEN '
badquerycand += str(ra1)
badquerycand += ' AND '
badquerycand += str(ra2)
badquerycand += ' AND c.decl BETWEEN '
badquerycand += str(decl1)
badquerycand += ' AND '
badquerycand += str(decl2)
badquerycand += ' AND o.jdmin <'
badquerycand += str(startdate)
badquerycand += ' AND '
badquerycand += candquality
#SQL query
print(badquerycand)

#Adds all bad candidate dictionaries to a list
bad_candid_results = mycursor.execute(badquerycand)
n = mycursor.rowcount
print('found %d potential candidates from bad objects' % n)
bad_candids = mycursor.fetchall()

print("-----------------------------------------------------------------------------------")

"""Finds all transients within ra, dec limits where transient jd is after cutoff"""
#all candidate query
querycand = 'SELECT candid, jd, ra, decl, magpsf, fid FROM candidates c WHERE '
querycand += ' c.ra BETWEEN '
querycand += str(ra1)
querycand += ' AND '
querycand += str(ra2)
querycand += ' AND c.decl BETWEEN '
querycand += str(decl1)
querycand += ' AND '
querycand += str(decl2)
querycand += ' AND c.jd BETWEEN '
querycand += str(startdate)
querycand += ' AND '
querycand += str(enddate)
querycand += ' AND '
querycand += candquality

#SQL query
print(querycand)

#All candidate dictionaries added to list
all_candid_results = mycursor.execute(querycand)
k = mycursor.rowcount
print('found %d candidates' % k)
all_candids = mycursor.fetchall()

print("-----------------------------------------------------------------------------------")

"""Finds all objects within ra, dec limits where jdmin is after cut-off."""
#object query
queryobj = 'SELECT objectId, ramean, decmean, rastd, decstd, maggmin, maggmean, maggmax, magrmin, magrmean, magrmax FROM objects WHERE '
queryobj += 'ramean BETWEEN '
queryobj += str(ra1)
queryobj += ' AND '
queryobj += str(ra2)
queryobj += ' AND decmean BETWEEN '
queryobj += str(decl1)
queryobj += ' AND '
queryobj += str(decl2)
queryobj += ' AND jdmin BETWEEN '
queryobj += str(startdate)
queryobj += ' AND '
queryobj += str(enddate)

#SQL query
print(queryobj)

#All object dictionaries added to list
obj_results = mycursor.execute(queryobj)
m = mycursor.rowcount
print('found %d objects' % m)
objects = mycursor.fetchall()

print("-----------------------------------------------------------------------------------")

mycursor.close()


# In[5]:


"""Removes any entries in all_candidates list which also occur in the bad_candidates list. 
This means that any transients observed are actually due to an object who's earliest detection was before the cut-off
and so cannot be considered."""
expired_candids = []
for i in all_candids: 
    if i in bad_candids:
        expired_candids.append(i)
good_candids = [i for i in all_candids if i not in expired_candids]

# should be less than all_candidates query
print(len(good_candids))


# In[6]:


"""Re-formats data to suitable form for matplotlib"""

#empty lists used for plotting
candid_list = []
cand_ra_list = []
cand_decl_list = []
cand_magpsf = []
cand_fid = []
cand_jd = []

obj_list = []
obj_ramean_list = []
obj_decmean_list = []
obj_rastd_list = []
obj_decstd_list = []
obj_maggmin_list = []
obj_maggmean_list = []
obj_maggmax_list = []
obj_magrmin_list = []
obj_magrmean_list = []
obj_magrmax_list = []


#separates data stored in dictionaries into various lists/numpy arrays 
for i in range(len(good_candids)):
    candid_list.append(good_candids[i]['candid'])
    cand_ra_list.append(good_candids[i]['ra'])
    cand_decl_list.append(good_candids[i]['decl'])
    cand_magpsf.append(good_candids[i]['magpsf'])
    cand_fid.append(good_candids[i]['fid'])
    cand_jd.append(good_candids[i]['jd'])
    
for j in range(len(objects)):
    obj_list.append(objects[j]['objectId'])
    obj_ramean_list.append(objects[j]['ramean'])
    obj_decmean_list.append(objects[j]['decmean'])
    obj_rastd_list.append(objects[j]['rastd'])
    obj_decstd_list.append(objects[j]['decstd'])
    obj_maggmin_list.append(objects[j]['maggmin'])
    obj_maggmean_list.append(objects[j]['maggmean'])
    obj_maggmax_list.append(objects[j]['maggmax'])
    obj_magrmin_list.append(objects[j]['magrmin'])
    obj_magrmean_list.append(objects[j]['magrmean'])
    obj_magrmax_list.append(objects[j]['magrmax'])

print("Done.")


# In[7]:


"""Finds associated probability for each candidate found within selected region"""

# lists for table plotting
cand_theta = []
cand_phi = []
cand_pix = []
cand_prob = []
cand_dist = []
cand_sigmadist = []
mag1 = []
mag2 = []

# finds theta and phi of each candidate
for i in cand_decl_list:
    cand_theta.append(np.radians(-i+90))
for i in cand_ra_list:
    cand_phi.append(np.radians(i))

# finds pixel number for each candidate
for i in range(len(candid_list)):
    cand_pix.append(hp.pixelfunc.ang2pix(NSIDE, cand_theta[i], cand_phi[i]))

# finds probability, distance and distance error for each candidate
for i in cand_pix:
    cand_prob.append(hpdata[i])
    cand_dist.append(distmu[i])
    cand_sigmadist.append(distsigma[i])

# separates each candidate's magnitude by r or g filters and appending to the two different lists.
for i in range(len(cand_pix)):
    if cand_fid[i] == 1:
        mag1.append(cand_magpsf[i])
        mag2.append(np.nan)
    if cand_fid[i] == 2:
        mag1.append(np.nan)
        mag2.append(cand_magpsf[i])
    
# finds the modified jd date of candidate
for i in range(len(cand_jd)):
    cand_jd[i] = cand_jd[i] - 2400000.5


# In[8]:


# lists for table plotting
obj_meantheta = []
obj_thetaerr = []
obj_meanphi = []
obj_phierr = []
obj_prob = []
obj_dist = []
obj_sigmadist = []
obj_pix = []

# finds theta and phi, means and errors of each object
for i in obj_decmean_list:
    obj_meantheta.append(np.radians(-i+90))
    obj_thetaerr.append(np.radians(i))
for i in obj_ramean_list:
    obj_meanphi.append(np.radians(i))
    obj_phierr.append(np.radians(i))

# determines object's pixel index
for i in range(len(obj_list)):
    obj_pix.append(hp.pixelfunc.ang2pix(NSIDE, obj_meantheta[i], obj_meanphi[i]))
    
# determines object's probability, distance and distance error
for i in obj_pix:
    obj_prob.append(hpdata[i])
    obj_dist.append(distmu[i])
    obj_sigmadist.append(distsigma[i])


# In[9]:


merg_dist = 39.5 #known distance to merger in Mpc
t = Table.read('refined_gw170817_phot_data.csv', format = 'ascii.csv')

derived_abs_mag = []
for i in t['magnitude']:
    derived_abs_mag.append(i-5*np.log10(merg_dist*1000000/10))
derived_abs_magmean = np.mean(derived_abs_mag)


# In[10]:


derived_app_mag1 = []
derived_app_mag2 = []
for i in range(len(cand_pix)):
    if cand_fid[i] == 1:
        derived_app_mag1.append(derived_abs_magmean + 5*np.log10(mag1[i]*1000000/10))
        derived_app_mag2.append(np.nan)
        
    if cand_fid[i] == 2:
        derived_app_mag2.append(derived_abs_magmean + 5*np.log10(mag2[i]*1000000/10))
        derived_app_mag1.append(np.nan)


# In[11]:


mag_diff = []
for i in range(len(candid_list)):
    if np.isnan(mag1[i]) == True:
        mag_diff.append(derived_app_mag2[i] - mag2[i])
    else:
        mag_diff.append(derived_app_mag1[i] - mag1[i])


# In[22]:


# identifies suitable columns for table
name = candid_list
mjd = cand_jd
ra = cand_ra_list
dec = cand_decl_list
magg = mag1
magr = mag2
prob = cand_prob
distance = cand_dist
sigmadistance = cand_sigmadist
deriv_appgmag = derived_app_mag1
deriv_apprmag = derived_app_mag2
# creates table
t1 = Table([name, mjd, ra, dec, magg, magr, distance, sigmadistance, prob, deriv_appgmag, deriv_apprmag, mag_diff], names=('Candidate number', 'MJD', 'RA', 'Dec', 'g-Magnitude', 'r-Magnitude', 'Distance, d (Mpc)', 'sigmad', 'Probability', 'Derived g-Magnitude', 'Derived r-Magnitude', 'mag diff'), meta={'name': 'candidate table'})

#order's table by highest probability
cand_idx_prob = np.where(t1['Probability'] >= max(hpdata)/1000)
candid_table = t1[cand_idx_prob]
candid_table.sort('Probability')

candid_table.reverse()

# shows table
candid_table.show_in_notebook()



# In[13]:


# identifies suitable columns for table
name = obj_list
ra = obj_ramean_list
rastd = obj_rastd_list
dec = obj_decmean_list
decstd = obj_decstd_list
maggmin = obj_maggmin_list
maggmean = obj_maggmean_list
maggmax = obj_maggmax_list
magrmin = obj_magrmin_list
magrmean = obj_magrmean_list
magrmax = obj_magrmax_list
prob = obj_prob
distance = obj_dist
sigmadistance = obj_sigmadist

# creates table
t2 = Table([name, ra, rastd, dec, decstd, maggmin, maggmean, maggmax, magrmin, magrmean, magrmax, distance, sigmadistance, prob], names=('Object ID', 'mean RA', 'RA error', 'mean Dec', 'Dec error', 'min g-magnitude', 'mean g-magnitude', 'max g-magnitude', 'min r-magnitude', 'mean r-magnitude', 'max r-magnitude', 'Distance, d (Mpc)', 'sigmad', 'Probability'), meta={'name': 'candidate table'})
obj_idx_prob = np.where(t2['Probability'] >= max(hpdata)/1000)
obj_table = t2[obj_idx_prob]

#order's table by highest probability
obj_table.sort('Probability')
obj_table.reverse()

# shows table
obj_table.show_in_notebook()


# In[15]:


"""plots separated data on an RA VS Dec scatter graph. Blue crosses represent a candidate, 
black spots represent an object (with error)"""

#scatters candidates
plt.scatter(cand_ra_list,cand_decl_list, s=1, marker = 'x')

#scatters objects
###for k in range(len(obj_list)):
###    plt.errorbar(obj_ramean_list[k], obj_decmean_list[k], yerr=obj_rastd_list[k], xerr=obj_decstd_list[k], ecolor = 'r', elinewidth = 0.5, capsize = 0.5)
###    plt.scatter(obj_ramean_list[k], obj_decmean_list[k], s = 5, color = 'k')
    
#graph formatting
plt.title('Candidates within GW170817')
plt.xlabel('RA')
plt.ylabel('Declination')

plt.show

