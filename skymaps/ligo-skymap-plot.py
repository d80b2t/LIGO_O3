



from ligo.skymap.tool.ligo_skymap_plot import main
from astropy.utils.data import download_file

#filename = download_file('https://dcc.ligo.org/public/0122/P1500227/012/bayestar_gstlal_C01.fits.gz', cache=True)

main([filename, '--annotate', '--contour', '50', '90'])
ligo-skymap-plot bayestar_gstlal_C01.fits.gz -o bayestar.png --annotate --contour 50 90
