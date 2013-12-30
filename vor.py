'''
(c) Copyright 2014 William Seviour. All Rights Reserved.
     
This file is part of vortex-moments. Please see README.md for more
information, including citations. 
 
vortex-moments is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
 
vortex-moments is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.
 
You should have received a copy of the GNU General Public License
along with vortex-moments.  If not, see <http://www.gnu.org/licenses/>.
'''

import numpy as np
from scipy.interpolate import griddata  

#################### GLOBAL PARAMETERS ##################
A = 6374.0e3       # Earth radius
DEGRAD = np.pi/180. # conversion of degrees to radians
RADDEG = 180./np.pi # conversion of radians to degrees 


   
def calc_moments(field, lats, lons, hemisphere='NH', field_type='GPH', \
                 edge=3.02e4):
    """Main procedure for calculating vortex moments
    
    **Arguments**
    
    *field*
        A :py:class:`numpy.ndarray` or :py:class:`numpy.ma.core.MasekdArray`
        with two dimensions. The 0th dimension should represent latitude, and 
        the 1st dimension longitude.
        
     *lats*
        A one dimensional array or list with the latitude grid values for field
    
     *lons*
        A one dimensional array or list with the longitude grid values for field
        
     **Optional Arguments:**
     
     *hemisphere* 
        A string either 'NH' or 'SH' which sets the calculation to the northern 
        or southern hemisphere respectively. Defaults to northern hemisphere.
        
     *field type*
        A string either 'GPH' or 'PV' which allows moments to be calculated over 
        geopotential heigh or potential vorticity respectively. Defaults to 
        geopotential height
        
     *edge*
        Value at the vortex edge. Defaults to 3.02e4m which is the climatological 
        value of zonal-mean geopotential height at 10hPa, 60N in ERA-Interim.
        
     **Returns:**
        A dictionary of moment diagnostic values.  
            
    """
    field_hem, lats_hem = select_hemisphere(field,lats,hemisphere)
    field_cart, x, y = sph_to_car(field_hem,lons,lats_hem,hemisphere)
    field_vtx = isolate_vortex(field_cart, edge, field_type)
    angle, aspect_ratio, objective_area, kurtosis, latcent, loncent = \
            moment_integrate(field_vtx, x, y,edge)
        
    return {'angle':angle, 'aspect_ratio':aspect_ratio, \
            'objective_area':objective_area, 'kurtosis':kurtosis, \
            'centroid_latitude':latcent, 'centroid_longitude':loncent}



def select_hemisphere(field, lats, hemisphere):
    """
    Selects appropriate hemisphere in field and lats
    """
    
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

    for ilon in range(len(lons)): # -1s needed?
        for ilat in range(len(lats)):
            
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
                    xypoints = np.vstack((xypoints,np.array((x[ilon,ilat], \
                                                             y[ilon,ilat]))))
                 
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
                    xypoints = np.vstack((xypoints,np.array((x[ilon,ilat], \
                                                             y[ilon,ilat]))))
                                  
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
    """
    Performs moment diagnostic calculations on cartesian field
    """

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
    loncent = np.arctan(centx/centy)*RADDEG
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

    # Calculate angle between x-axis and major axis of ellipse                  
    angle = 0.5*np.arctan((2*J11)/(J20-J02))*RADDEG

    # Calculate aspect ratio
    aspect_ratio = np.sqrt(abs(( (J20+J02) + np.sqrt(4*(J11**2)+(J20-J02)**2) ) / \
                               ( (J20+J02) - np.sqrt(4*(J11**2)+(J20-J02)**2) ))) 
    ar = aspect_ratio #short name for later calculation

    # Calculate objective area
    objective_area = Marea/edge 


    # Calculate excess kurtosis 
    kurtosis = M00 * (J40+2*J22+J04)/((J20+J02)**2) \
                         - (2./3.)*( (3*(ar**4)+2*(ar**2)+3) / (((ar**2)+1)**2) ) 
                   
    return angle, aspect_ratio, objective_area, kurtosis, latcent, loncent 
               
    
    
    
    
                  
                
                
