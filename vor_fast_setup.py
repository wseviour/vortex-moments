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
    nlon_nlat = lons.shape[0]*lats.shape[0]
    
    if hemisphere == 'NH':
        x = (np.cos(lons[:,np.newaxis])*np.cos(lats))/ (1. + np.sin(lats[np.newaxis,:]))               
        y = (np.sin(lons[:,np.newaxis])*np.cos(lats[np.newaxis,:]))/(1. + np.sin(lats[np.newaxis,:]))
          
    elif hemisphere == 'SH':
        x = (np.cos(lons[:,np.newaxis])*np.cos(lats[np.newaxis,:]))/(1. - np.sin(lats[np.newaxis,:]))               
        y = (-np.sin(lons[:,np.newaxis])*np.cos(lats[np.newaxis,:]))/(1. - np.sin(lats[np.newaxis,:]))     
        
    else:
        raise ValueError()        
    
    xypoints = np.stack([x.reshape(nlon_nlat),y.reshape(nlon_nlat)], axis = 1)

    return xypoints
