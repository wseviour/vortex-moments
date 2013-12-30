"""
This example uses vor.py to plot geopotential height derived moment diagnostics 
for Feb-Mar 1979. This covers the splitting event around 21st February

** Using a faster method ** (takes approx. 10 mins)
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import vor_fast
import vor_fast_setup

# Read in NetCDF file with geopotential height values
ncin = Dataset('examples_data/GPH_10hPa_1979_FebMar_daily.nc', 'r')
gph = ncin.variables['gph'][:].squeeze()
lons = ncin.variables['longitude'][:]
lats = ncin.variables['latitude'][:]
days = ncin.variables['time'][:]
ncin.close()

# Set up cartesian mapping xypoints and restrict to NH
gph_nh, lats_nh, xypoints = vor_fast_setup.setup(gph,lats,lons,'NH')

# Set up moment diagnostics
aspect = np.empty(0)
latcent = np.empty(0)


# Calculate diagnostics for each day
for iday in range(len(days)):
    print 'Calculating moments for day '+str(iday)
    moments = vor_fast.calc_moments(gph_nh[iday,:,:],lats_nh,lons,xypoints,'NH','GPH',3.02e4) 
    aspect = np.append(aspect, moments['aspect_ratio'])
    latcent = np.append(latcent, moments['centroid_latitude'])

# Plot timeseries    
plt.subplot(2,1,1)
plt.plot(aspect)
plt.title('Aspect ratio')
plt.subplot(2,1,2)
plt.plot(latcent) 
plt.title('Centroid latitude')
plt.show()
