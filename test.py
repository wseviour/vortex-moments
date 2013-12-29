import cdms2
import numpy as np
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import interp 
#http://matplotlib.org/basemap/users/mapcoords.html
#http://earthpy.org/interpolation_between_grids_with_basemap.html

f = cdms2.open('/home/will/data/GPH_10hPa_2000_daily_noleaps.nc', 'r')
gph = f['gph']

field = gph[0,0,0:80,:]

hemisphere = 'NH'
lons = field.getLongitude()[:] 
lats = field.getLatitude()[:] 

# are lons -180:180 or 0:360 ?
if max(lons) > 180:
    lon1 = lons.copy()
    for n, l in enumerate(lon1):
        if l >= 180:
            lon1[n]=lon1[n]-360. 
    lons = lon1

    half_lons = lons.shape[0]/2 

    field1 = field[:,0:half_lons]
    field2 = field[:,half_lons:]
    lon1 = lons[0:half_lons]
    lon2 = lons[half_lons:]
    field_new = np.hstack((field2, field1))
    lons_new = np.hstack((lon2, lon1))

# lats must be increasing
if lats[0] > lats[1]:
    lats_new = np.flipud(lats)
    field_new = np.flipud(field_new)
    
m = Basemap(projection='npstere', boundinglat=0, lon_0=0)
x, y = m(*np.meshgrid(lons_new,lats_new)) 

result = interp(field_new, lons_new, lats_new, y, x)#, checkbounds=False, masked=False, order=1)

