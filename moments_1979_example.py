"""
This example uses vor.py to plot geopotential height derived moment diagnostics 
for Feb-Mar 1979. This covers the splitting event around 21st February

WARNING: This takes about 15 mins to run
"""

from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import vor

# Read in NetCDF file with geopotential height values
ncin = Dataset('examples_data/GPH_10hPa_1979_FebMar_daily.nc', 'r')
gph = ncin.variables['gph'][:].squeeze()
lons = ncin.variables['longitude'][:]
lats = ncin.variables['latitude'][:]
days = ncin.variables['time'][:]
ncin.close()

# Set up moment diagnostics
aspect = np.empty(0)
latcent = np.empty(0)
kurtosis = np.empty(0)
obj_area = np.empty(0)

# Calculate diagnostics for each day
for iday in range(len(days)):
    print('Calculating moments for day '+str(iday))
    moments = vor.calc_moments(gph[iday,:,:],lats,lons,'NH','GPH',3.02e4) 
    aspect = np.append(aspect, moments['aspect_ratio'])
    latcent = np.append(latcent, moments['centroid_latitude'])
    kurtosis = np.append(kurtosis, moments['kurtosis'])
    obj_area = np.append(obj_area, moments['objective_area'])
    
# Plot timeseries    
fig = plt.figure(figsize = (12,9))
plt.subplot(4,1,1)
plt.plot(aspect)
plt.title('Aspect ratio')
plt.subplot(4,1,2)
plt.plot(latcent) 
plt.title('Centroid latitude')
plt.subplot(4,1,3)
plt.plot(kurtosis)
plt.title('Kurtosis')
plt.subplot(4,1,4)
plt.plot(obj_area)
plt.title('Objective area')
plt.xlabel('days from 1st Feb 1979')
plt.tight_layout()
plt.savefig('test.png')
