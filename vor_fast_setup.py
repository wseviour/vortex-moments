import numpy as np

#################### GLOBAL PARAMETERS ##################
A = 6374.0e3       # Earth radius
DEGRAD = np.pi/180. # conversion of degrees to radians
RADDEG = 180./np.pi # conversion of radians to degrees 

def setup(field,lats,lons,hemisphere='NH'):

    fieldh, latsh = select_hemisphere(field, lats, hemisphere)
    
    xypoints = calc_cart_map(lons,latsh,hemisphere)
    
    return fieldh, latsh, xypoints


def select_hemisphere(field, lats, hemisphere):
    """
    Selects appropriate hemisphere in field and lats
    """
    
    if hemisphere == 'NH':
        latsh = lats[np.where(lats > 20)]
        fieldh = field[:,np.where(lats > 20),:].squeeze()
    if hemisphere == 'SH':
        latsh = lats[np.where(lats < -20)]
        fieldh = field[:,np.where(lats < -20),:].squeeze()  
     
    return fieldh, latsh

# Note only data for correct hemisphere is used 
def calc_cart_map(lons, lats, hemisphere='NH'):
    # Convert lats and lons to radians
    lons = lons * DEGRAD
    lats = lats * DEGRAD
    x = np.zeros((len(lons), len(lats)))
    y = np.zeros((len(lons), len(lats)))
    try: 
        del xypoints
    except NameError:
        pass 


    for ilon in range(len(lons)):
        for ilat in range(len(lats)):
            
            if hemisphere == 'NH':
                x[ilon,ilat] = (np.cos(lons[ilon])*np.cos(lats[ilat]))/ \
                                   (1. + np.sin(lats[ilat]))               
                y[ilon,ilat] = (np.sin(lons[ilon])*np.cos(lats[ilat]))/ \
                                   (1. + np.sin(lats[ilat]))
                try:
                    xypoints
                except NameError:                   
                    xypoints = np.array((x[ilon,ilat],y[ilon,ilat]))
                else:
                    xypoints = np.vstack((xypoints,np.array((x[ilon,ilat], \
                                                             y[ilon,ilat]))))
                 
            if hemisphere == 'SH':
                x[ilon,ilat] = (np.cos(lons[ilon])*np.cos(lats[ilat]))/ \
                                   (1. - np.sin(lats[ilat]))               
                y[ilon,ilat] = (-np.sin(lons[ilon])*np.cos(lats[ilat]))/ \
                                   (1. - np.sin(lats[ilat]))            
                try:
                    xypoints
                except NameError:                   
                    xypoints = np.array((x[ilon,ilat],y[ilon,ilat]))
                else:
                    xypoints = np.vstack((xypoints,np.array((x[ilon,ilat], \
                                                             y[ilon,ilat]))))
    return xypoints
