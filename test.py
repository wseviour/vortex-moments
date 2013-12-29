import cdms2
import numpy as np
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import interp 
import matplotlib.pyplot as plt
from vor import sph_to_car
#http://matplotlib.org/basemap/users/mapcoords.html
#http://earthpy.org/interpolation_between_grids_with_basemap.html

f = cdms2.open('/home/will/data/GPH_10hPa_2000_daily_noleaps.nc', 'r')
gph = f['gph']

hemisphere = 'NH'
lons = gph.getLongitude()[:] 
lats = gph.getLatitude()[:] 
a = gph[0,0,:,:].data   

x = vor.calc_moments(a,lats,lons,'NH','GPH',29600) 


#field_car, x, y = sph_to_car(field, lons, lats, hemisphere)

#plt.subplot(121)
#plt.imshow(field_car)
#plt.subplot(122)
#plt.imshow(field)
#plt.show()

#from vor import isolate_vortex

#vtx_field = isolate_vortex(field_car, 29600, 'GPH')

#plt.imshow(vtx_field)    
#plt.show()

#from vor import calc_moments

#moments = calc_moments(vtx_field, x, y, 29600) 


