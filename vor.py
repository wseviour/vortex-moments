'''
HEADER INFO HERE

Inputs: (field[lat, lon], lats, lons, edgevalue )

'''
import numpy as np
from scipy.interpolate import griddata  

#################### GLOBAL PARAMETERS ##################
a = 6374.0e3       # Earth radius
degrad = np.pi/180 # conversion of degrees to radians




# Note only data for correct hemisphere is used 
def sph_to_car(field, lons, lats, hemisphere='NH'):
    """
    Convert field on sperical (lon,lat) coordinates cartesian coordinates.
    Uses a polar stereographic projection, where the hemisphere needs to be 
    specified. 
    
    """
    x = np.zeros((len(lons), len(lats)))
    y = np.zeros((len(lons), len(lats)))
    try: 
        del xypoints
    except NameError:
        pass 

    xyvals = np.empty(0)

    for ilon in range(len(lons)-1): # -1s needed?
        for ilat in range(len(lats)-1):
            
            if hemisphere == 'NH':
                x[ilon,ilat] = (np.cos(lons[ilon])*np.cos(lats[ilat]))/ \
                                   (1. + np.sin(lats[ilat]))               
                y[ilon,ilat] = (np.sin(lons[ilon])*np.cos(lats[ilat]))/ \
                                   (1. + np.sin(lats[ilat]))
                xyvals = np.append(xyvals, field[ilat,ilon])
                try:
                    xypoints
                except NameError:                   
                    xypoints = np.array((x[ilon,ilat],y[ilon,ilat]))
                else:
                    xypoints = np.vstack((xypoints,np.array((x[ilon,ilat],y[ilon,ilat]))))
                 
            if hemisphere == 'SH':
                x[ilon,ilat] = (np.cos(lons[ilon])*np.cos(lats[ilat]))/ \
                                   (1. - np.sin(lats[ilat]))               
                y[ilon,ilat] = (-np.sin(lons[ilon])*np.cos(lats[ilat]))/ \
                                   (1. - np.sin(lats[ilat]))            
                xyvals = np.append(xyvals, field[ilat,ilon])
                try:
                    xypoints
                except NameError:                   
                    xypoints = np.array((x[ilon,ilat],y[ilon,ilat]))
                else:
                    xypoints = np.vstack((xypoints,np.array((x[ilon,ilat],y[ilon,ilat]))))
                                  
    cart_x_points = -1.+np.arange(len(lons))/(0.5*len(lons))             
    cart_y_points = -1.+np.arange(len(lons))/(0.5*len(lons))
    cart_gridx, cart_gridy = np.meshgrid(cart_x_points,cart_y_points)

    field_cart = griddata(xypoints, xyvals, (cart_gridx,cart_gridy), \
                        method='linear')  # Might want to change to cubic etc?
                              
    return field_cart
                    
    
                
                
                
                
                
                
                
