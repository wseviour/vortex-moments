'''  
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


   
def calc_moments(field, lats, lons, xypoints, hemisphere='NH', field_type='GPH', \
                 edge=3.02e4, resolution='full'):
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
     
     *xypoints*
        Mapping of lat/lon points to cartesian as calculated by vor_fast_setup.py
        
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
    print 'Calculating for resolution: '+resolution
    field_cart, x, y = sph_to_car(field,lons,lats,xypoints,resolution)
    field_vtx = isolate_vortex(field_cart, edge, field_type)
   
    aspect_ratio, latcent = moment_integrate(field_vtx, x, y,edge)
        
    return {'aspect_ratio':aspect_ratio, 'centroid_latitude':latcent}

                                  
def sph_to_car(field, lons, lats, xypoints,resolution):
    nlons = len(lons)
    nlats = len(lats)
    xyvals = []
    for ilon in range(nlons): # -1s needed?
        for ilat in range(nlats):
            xyvals.append(field[ilat,ilon])

    if resolution == 'full':       
        cart_x_points = -1.+np.arange(nlons)/(0.5*nlons)             
        cart_y_points = -1.+np.arange(nlons)/(0.5*nlons)
    elif resolution == 'low':
        cart_x_points = -1.+np.arange(50)/(0.5*50)
        cart_y_points = -1.+np.arange(50)/(0.5*50)
    else:
        raise ValueError()
    
    cart_gridx, cart_gridy = np.meshgrid(cart_x_points,cart_y_points)

    field_cart = griddata(xypoints, np.array(xyvals), (cart_gridx,cart_gridy), \
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

    # Set up coords
    x = x[:,np.newaxis]
    y = y[np.newaxis,:]
    # Integrate over vortex
    M00 = np.sum(np.abs(vtx_field-edge)*(x**0)*(y**0))
    M10 = np.sum(np.abs(vtx_field-edge)*(x**1)*(y**0))
    M01 = np.sum(np.abs(vtx_field-edge)*(x**0)*(y**1))
    Marea = np.sum(np.abs(vtx_field-edge)*(x**0)*(y**0)*box_area)
            
    # Calculate centroid 
    centx = M10/M00
    centy = M01/M00

    # Convert back to polar coordinates 
    R = centx**2 + centy**2
    #centroid latitude
    latcent = np.arcsin((1-R)/(1+R))*RADDEG

    #Intergrate relative moment diagnostics 
    J11 = np.sum(np.abs(vtx_field-edge)*((x-centx)**1)*((y-centy)**1))
    J20 = np.sum(np.abs(vtx_field-edge)*((x-centx)**2)*((y-centy)**0))
    J02 = np.sum(np.abs(vtx_field-edge)*((x-centx)**0)*((y-centy)**2)) 


    # Calculate aspect ratio
    aspect_ratio = np.sqrt(abs(( (J20+J02) + np.sqrt(4*(J11**2)+(J20-J02)**2) ) / \
                               ( (J20+J02) - np.sqrt(4*(J11**2)+(J20-J02)**2) ))) 

                   
    return aspect_ratio, latcent 
               
    
    
    
    
                  
                
                
