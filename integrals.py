# -*- coding: utf-8 -*-
"""
Integrals of the vertex and segment distance functions required for
calculating the average distance from one polyline to another.

Created on Fri Oct 27 14:56:55 2023
@author: Barry Kronenfeld
"""
import math
import utils_hausdorff

def component_integral(L,k0,k1,dr):
    if dr[0]: # component is a segment
        return segment_integral(L,k0,k1,dr)
    else:
        return vertex_integral(L,k0,k1,dr)
            

def vertex_integral(L,k0,k1,dr):
    """
        Calculates the area under the vertex distance function to a vertex on
        polyline B between two points along a segment on polyline B
        Parameters
        ----------
        L : float
            length of segment
        k0 : float
            proportional distance along segment to first point
        k1 : float 
            proportional distance along segment two second point
        dist_rep : distance representation
            distance representation of vertex
        Returns
        ----------
        float
            The area under the distance function from the segment to the vertex
            between the two points.

    """
    kv = dr[1]
    q = dr[2]
    start_offset = k0 - kv
    finish_offset = k1 - kv
    start_dist = ((start_offset)**2 + q**2)**0.5    

    if q == 0 or start_dist + start_offset == 0:
        # point is on line segment within floating point precision error
        # calculate based on straight line
        area = (L**2) * ((start_offset**2) + (finish_offset**2)) / 2
    else:
        finish_dist = ((finish_offset)**2 + q**2)**0.5    
        term1 = finish_offset * finish_dist
        term2 = start_offset * start_dist        
        term3 = (q**2) * math.log((finish_dist + finish_offset)/(start_dist + start_offset))
        area = (L**2) * (term1 - term2 + term3) / 2
    return area
    
def segment_integral(L,k0,k1,dr): 
    """
        Calculates the area under the segment distance function to a segment on
        polyline B between two points along a segment of polyline A
        
        Parameters
        ----------
        L : float
            length of segment on polyline A
        k0 : float
            proportional distance along segment to first point
        k1 : float 
            proportional distance along segment two second point
        ks : float
            proportional distance along segment to intersection with segment of 
            polyline B, or None if two segments are parallel
        cos_theta : float
            cosine of angle between two segments, or q-distance (proportional 
            distance) between segments if they are parallel.
        Returns
        ----------
        float
            The area under the distance function 
    """
    ks = dr[1]
    cos_theta = dr[2]
    # handle special cases
    if ks == None:
        return abs(cos_theta * (k1 - k0) * L**2)
    else:
        dist_rep = (True,ks,cos_theta)
        d0 = utils_hausdorff.segDistance(dist_rep, k0, L)
        d1 = utils_hausdorff.segDistance(dist_rep, k1, L)
        if (k0 < ks) == (k1 < ks):
            return abs(0.5 * (d0 + d1) * (k1 - k0) * L)
        else:
            return abs(0.5 * d0 * (k0 - ks) * L) + abs(0.5 * d1 * (k1 - ks) * L)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    