'''
HEADER INFO HERE

Inputs: (field[lat, lon], lats, lons, edgevalue )

'''
import numpy as np
from scipy.interpolate import griddata  

#################### GLOBAL PARAMETERS ##################
A = 6374.0e3       # Earth radius
DEGRAD = np.pi/180. # conversion of degrees to radians
RADDEG = 180./np.pi # conversion of radians to degrees 


def calc_moments(field, lats, lons, hemisphere, field_type, edge):

    field_hem, lats_hem = select_hemisphere(field,lats,hemisphere)
    field_cart, x, y = sph_to_car(field_hem,lons,lats_hem,hemisphere)
    field_vtx = isolate_vortex(field_cart, edge, field_type)
    angle, aspect_ratio, equivalent_area, kurtosis, latcent, loncent = \
            moment_integrate(field_vtx, x, y,edge)
        
    return angle, aspect_ratio, equivalent_area, kurtosis, latcent, loncent



def select_hemisphere(field, lats, hemisphere):
    if hemisphere == 'NH':
        latsh = lats[np.where(lats > 0)]
        fieldh = field[np.where(lats > 0),:].squeeze()
    if hemisphere == 'SH':
        latsh = lats[np.where(lats < 0)]
        fieldh = field[np.where(lats < 0),:].squeeze()  
     
    return fieldh, latsh

# Note only data for correct hemisphere is used 
def sph_to_car(field, lons, lats, hemisphere='NH'):
    """
    Convert field on sperical (lon,lat) coordinates cartesian coordinates.
    Uses a polar stereographic projection, where the hemisphere needs to be 
    specified. 
    
    """
    # Convert lats and lons to radians
    lons = lons * DEGRAD
    lats = lats * DEGRAD
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
                              
    return field_cart, cart_x_points, cart_y_points


 
    
    
    
def isolate_vortex(field_cart, edge, field_type='GPH'):
    """
    Replace vortex into region (outside) equal to field on vortex edge, and 
    region (inside) with values of vortex. 
    For GPH vortex is less than surrounding, for PV it is greater. 
    """
    if field_type == 'GPH':
        field_cart[np.where(field_cart > edge)] = edge
        field_cart[np.isnan(field_cart)] = edge # set NaN regions to edge 
    elif field_type == 'PV':
        field_cart[np.where(field_cart < edge)] = edge
        field_cart[np.isnan(field_cart)] = edge # set NaN regions to edge
    else:
        raise ValueError() # maybe want more specific error type

    return field_cart




def moment_integrate(vtx_field, x, y,edge):
    # x and y are cartesian gridpoints; vtx_field is cartesian field
    # edge is value on vortex edge 

    box_length = 2*A/len(x)
    box_area = box_length**2

    # Set up moment diagnostics
    M00 = 0
    M10 = 0
    M01 = 0
    Marea = 0
    # Integrate over vortex
    for ix in range(len(x)):
        for iy in range(len(y)):
            
            M00 += abs(vtx_field[ix,iy]-edge)*(x[ix]**0)*(y[iy]**0) 
            M10 += abs(vtx_field[ix,iy]-edge)*(x[ix]**1)*(y[iy]**0) 
            M01 += abs(vtx_field[ix,iy]-edge)*(x[ix]**0)*(y[iy]**1)
            Marea += abs(vtx_field[ix,iy]-edge)*(x[ix]**0)*(y[iy]**0)*box_area
            
    # Calculate centroid 
    centx = M10/M00
    centy = M01/M00

    # Convert back to polar coordinates 
    R = centx**2 + centy**2
    loncent = 90. - np.arctan(centx/centy)*RADDEG
    latcent = np.arcsin((1-R)/(1+R))*RADDEG

    # Set up relative moment diagnostics 
    J11=0
    J20=0
    J02=0
    J22=0
    J40=0
    J04=0
    for ix in range(len(x)):
        for iy in range(len(y)):
            
            J11 += abs(vtx_field[ix,iy]-edge)*((x[ix]-centx)**1)*((y[iy]-centy)**1)
            J20 += abs(vtx_field[ix,iy]-edge)*((x[ix]-centx)**2)*((y[iy]-centy)**0)
            J02 += abs(vtx_field[ix,iy]-edge)*((x[ix]-centx)**0)*((y[iy]-centy)**2)
            J22 += abs(vtx_field[ix,iy]-edge)*((x[ix]-centx)**2)*((y[iy]-centy)**2)
            J40 += abs(vtx_field[ix,iy]-edge)*((x[ix]-centx)**4)*((y[iy]-centy)**0)
            J04 += abs(vtx_field[ix,iy]-edge)*((x[ix]-centx)**0)*((y[iy]-centy)**4)   

    # Calculate angle between x-axis and majoe axis of ellipse
    angle = 0.5*np.arctan((2*J11)/(J20-J02))                

    # Calculate aspect ratio
    aspect_ratio = np.sqrt(abs(( (J20+J02) + np.sqrt(4*(J11**2)+(J20-J02)**2) ) / \
                               ( (J20+J02) - np.sqrt(4*(J11**2)+(J20-J02)**2) ))) 
    ar = aspect_ratio #short name for later calculation

    # Calculate equivalent area
    equivalent_area = Marea/edge 


    # Calculate excess kurtosis 
    kurtosis = M00 * (J40+2*J22+J04)/((J20+J02)**2) \
                         - (2./3.)*( (3*(ar**4)+2*(ar**2)+3) / (((ar**2)+1)**2) ) 

    moments = {'phi_angle':angle, 'aspect_ratio':aspect_ratio, \
               'equivalent_area':equivalent_area, 'kurtosis':kurtosis, \
               'centroid_latitude':latcent, 'centroid_longitude':loncent}
               
    return angle, aspect_ratio, equivalent_area, kurtosis, latcent, loncent 
               
    
    
    
    
                  
                
                
