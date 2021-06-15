# -*- coding: utf-8 -*-
"""
Calculates the Hausdorff distance between two polylines
following the method described in:
J.F. Hangouet (1995). Computation of the Hausdorff distance between plane 
vector polylines. Proceedings, Auto Carto 12, pp. 1-10.

A joint project of the Spring 2021 semester GEO 4910 class at Eastern Illinois 
University.

Team members:
    Luke Jansen
    Tanner Jones
    Farouk Olaitan
    Megshi Thakur
    Barry Kronenfeld (Instructor)
    
License: MIT License

If you publish or create derivate software using this code, acknowledgement 
would be appreciated.

********************************************************************

CONVENTIONS:
    The two polylines are referred to as A and B.

    Lowercase a or b indicates the index of either a vertex or segment on 
    A or B.

    For a polyline with n vertices, vertices are indexed 0 to n-1 and segments 
    are indexed 0 to n-2

DEFINITIONS:
    k-value:        float
        The "relative location" along a segment of polyline A.
        - At the start of the segment, k=0
        - At the end of the segment, k=1
    q-value:        float
        The "relative distance" to a segment of polyline A, i.e. the 
        perpendicular distance to the segment divided by the length of the 
        segment.
    component:      (bool, int)
        A tuple of two values identifying a component of polyline B, consisting of:
            bool:   If True the component is a segment, if False it is a vertex
            int:      The index of the vertex or segment on polyline B
    distance representation:      (bool, float, float)
        A tuple of three values containing information needed to compute the 
        distance from any point on a segment "a" of A to component b of B:
            bool:   If True b is a segment, if False it is a vertex
            float:  The k-value of the projection of vertex b onto segment a, or
                    The k-value of the point of intersection of segment b with segment a, or
                    None if b is a segment parallel to a
            float:  The sin of the angle between a & b (if b is a non-parallel segment)
                    The q-distance from b to a (if b is a vertex or segment parallel to a)
"""

import utils_geom as g
import utils_hausdorff as hu

def polyline_hausdorff(A,B):
    """
    Computes the Hausdorff distance between two polylines, i.e. the maximum
    distance from any point on either polyline to the nearest point on the
    other polyline.

    Parameters
    ----------
    A : [(x,y),...] list of tuples 
        The coordinates of a polyline.
    B : [(x,y),...] list of tuples 
        The coordinates of another polyline.
    
    Returns
    ----------
    h_dist : float
        The Hausdorff distance between the polylines.
    srcloc : (float,float)
        The location of the source of the Hausdorff distance.
        The source location is the location on either polyline that is furthest 
        from the nearest location on the other polyline.
    srcline : int
        Flag to indicate the source polyline. 
        If 0, the source location is on polyline A. 
        If 1, the source location is on polyline B.
    srccomp : (bool,int)
        The index of the component (vertex or segment) of the source polyline 
        containing the source location. 
        If bool is True then srccomp is a segment, otherwise it is a vertex.
    trgcomp1 : (bool, int)
        Target component #1.
        The index of the component of the other ("target") polyline that is 
        closest to the source location.
    trgcomp2 : (bool, int)
        Target component #2.
        A second component of the target polyline that is equally close to 
        the source location, or None if only one target component was found. 
        If the source component is a segment then then there will usually be 
        two equidistant nearest locations on the target polyline.    
    """
    # check in both directions
    h_dist_a,srcloca,srccompa,trgcompa1,trgcompa2 = hausdorff_unidirectional(A,B)
    h_dist_b,srclocb,srccompb,trgcompb1,trgcompb2 = hausdorff_unidirectional(B,A)
    # determine which is source and which is target    
    if h_dist_a > h_dist_b:
        srcloc = srcloca
        trgline = B
        trgcomp1 = trgcompa1
        trgcomp2 = trgcompa2        
    else:
        srcloc = srclocb
        trgline = A
        trgcomp1 = trgcompb1
        trgcomp2 = trgcompb2
    # determine target location(s)
    trglocs = []
    trglocs.append(hu.nearLoc(srcloc, trgline, trgcomp1))
    if trgcomp2 != None:        
        trglocs.append(hu.nearLoc(srcloc, trgline, trgcomp2))
    return h_dist_a,srcloc,trglocs


def hausdorff_unidirectional(A,B):
    """
    Computes the one-way Hausdorff distance from A to B, i.e. the maximum
    distance from any point on A to the nearest point on B.

    Parameters
    ----------
    A : [(x,y),...] list of tuples 
        The coordinates of the main polyline.
    B : [(x,y),...] list of tuples 
        The coordinates of the other polyline.
    Returns
    ----------
    H : float
        The unidirectional Hausdorff distance between the polylines.
    loc : (float, float)
        The location on A that is the source of the Hausdorff distance.
    srcComp : (bool,int)
        The component of A that is the source of the Hausdorff distance.
    trgComp1 : (bool, int)
        A component of B that is the target of the Hausdorff distance.
    trgComp2 : (bool, int)
        Another component of B that is the target of the Hausdorff distance, or
        None if there is only one target component.
    """
    # initialize component and distance lists
    vNearComp = []
    vertDist = []
    # get distance from each vertex on A to its nearest component on B
    for a in range(len(A)):
        b = hu.nearSegment(A, B, a)
        comp, d = hu.nearComponent(A, B, a, b)
        vNearComp.append(comp)
        vertDist.append(d)

    segsToCheck = [a for a in range(len(A)-1) if hu.checkSegment(vNearComp[a],vNearComp[a+1])]
    # get distance from segments on A to nearest components on B
    d_k_comps = []
    for a in segsToCheck:
        comp_list = hu.candidateComponents(A,B,a)
        d,k,comp1,comp2 = hausdorff_segment(A, B, a, vNearComp[a],vertDist[a],vNearComp[a+1],vertDist[a+1], comp_list)
        d_k_comps.append((d,k,comp1,comp2))
    # Hausdorff distance is maximum of distances from vertices and segments
    maxVertDist = max(vertDist)
    maxSegInfo = max(d_k_comps, key = lambda x: x[0])
    if maxVertDist > maxSegInfo[0]:
        a = vertDist.index(maxVertDist)
        loc = A[a]
        srcComp = (False,a)
        trgComp1 = vNearComp[a]
        trgComp2 = None
        H = maxVertDist
    else:
        a = segsToCheck[d_k_comps.index(maxSegInfo)]
        d,k,trgComp1,trgComp2 = maxSegInfo
        loc = g.location(A,a,k)
        srcComp = (True,a)
        H = d
    return H,loc,srcComp,trgComp1,trgComp2

def hausdorff_segment(A,B,a,startComp,startDist,endComp,endDist,comp_list,verbose = False,max_iterations = 1000000):
    """
    Computes the furthest distance from any point on segment a of polyline A 
    to the nearest point on polyline B

    Parameters
    ----------
    A : [(x,y),...] list of tuples 
        The coordinates of the main polyline.
    B : [(x,y),...] list of tuples 
        The coordinates of the other polyline.
    a : int
        The index of a segment on polyline A.
    startComp : comp
        Component of B nearest to vertex a.
    startDist : float
        Distance from vertex a to startComp
    endComp : comp
        Component of B nearest to vertex a+1.
    endDist : float
        Distance from vertex a+1 to endComp
    comp_list : [comp]
        List of all candidate components of B for the Hausdorff distance.
    verbose : bool
        If True, will print information on components traversed
    max_iterations : 1000000
        Flag to stop, in case of any unknown special cases that
        would cause an infinite loop.
    Returns
    ----------
    d : float
        The directed Hausdorff distance from segment a to polyline B.
    k : float
        The k-value of the source of the Hausdorff distance along segment a.        
    max_comp1 : (bool, int)
        First of (at least) two components on B that is the target of the Hausdorff distance
    max_comp2 : (bool, int)
        Second of (at least) two components on B that is the target of the Hausdorff distance
    """
    # We're going to "walk" along segment a from vertex a to vertex a+1,
    # keeping track of the nearest component on B as we go

    # make sure near component to start of A is in list
    if not startComp in comp_list:
        comp_list.append(startComp)
    
    # pre-calculate distance representations and effective intervals
    drs = [hu.distanceRepresentation(A, B, a, c) for c in comp_list]
    eis = [hu.effectiveInterval(A, B, a, c) for c in comp_list]
    
    # intitialize to start of segment a (k=0)
    prev_id = comp_list.index(startComp)
    ei = eis[prev_id]
    k = 0
    # in some cases, floating point precision errors will cause 
    # nearest component to have no effective interval.
    # for such cases, set the initial effective interval to [0,1]
    if ei==(float('-inf'),float('-inf')):
        ei = (0,1)
    # initialize maximum distance and related information to further of
    # start and end vertices of segment
    if startDist > endDist:
        max_comp1 = startComp    
        max_d = startDist
        max_k = 0 
    else:
        max_comp1 = endComp    
        max_d = endDist
        max_k = 1 
    max_comp2 = None
    
    if verbose:
        print("\nDistance to beat: {}".format(max_d))
        print("From k = {}".format(max_k))
        print("To {}".format(hu.component_label(max_comp1)))
    
    # keep going until end of line segment is reached
    numIts = 0
    while prev_id != -1 and numIts < max_iterations:
        numIts += 1
        # report component we're moving from
        if verbose:
            print("\nMoving from {}".format(hu.component_label(comp_list[prev_id])))
            print("incoming ei: {}".format(ei))
        # move to next component
        new_id,ei,d,k = _updateComponent(A,B,a,comp_list,eis,drs,prev_id,ei)
        # record components, distance at switch point
        if new_id != -1:
            if d > max_d:
                max_d = d
                max_k = k
                max_comp1 = comp_list[new_id]
                max_comp2 = comp_list[prev_id]
        # remove previous component from list
        del comp_list[prev_id]
        del drs[prev_id]
        del eis[prev_id]
        # move to next component
        if prev_id < new_id: # we deleted prev_id so subtract one from new_id
            prev_id = new_id - 1
        else:
            prev_id = new_id
        # report component we're moving to
        if verbose and new_id != -1:
            print("moving to {}".format(hu.component_label(comp_list[prev_id])))
            print("k: {}".format(k))
            print("ei: {}".format(ei))
            print("d: {}".format(d))

    if verbose:                
        # report final solution
        a1 = A[a]
        a2 = A[a+1]
        x = a1[0] + max_k * (a2[0]-a1[0])
        y = a1[1] + max_k * (a2[1]-a1[1])
        loc = (x,y)
        print("\nfurthest component(s) from {}: {},{}".format(hu.component_label((True,a)), hu.component_label(max_comp1), hu.component_label(max_comp2)))
        print("distance: {}".format(max_d))
        print("k-value: {}".format(max_k))
        print("location:\t {}\t{}".format(loc[0],loc[1]))
    return max_d, max_k, max_comp1, max_comp2

def _updateComponent(A,B,a,comps, eis, drs, prev_id, prev_ei):
    """
    Finds the next component of B while walking alon gsegment a of A. 
    Simultaneously updates the current effective interval

    Parameters
    ----------
    A :     [(x,y),...] list of tuples 
            The coordinates of the main polyline.
    B :     [(x,y),...] list of tuples 
            The coordinates of the other polyline.
    a :     int
            The index of a segment on polyline A.
    comps : [component]
            List of candidate components on B that could be target of Hausdorff distance.
    ei :    list of (float,float)
            K-values of effective intervals of above comps on segment a of A.
    drs :   list of (bool,float,float)
            Distance representations of above comps wrt seg a of A
    prev_id : int
            ID in comps of previous component of polyline B.
    prev_ei : (float,float)
            Previous effective interval on seg a before update

    Returns
    ----------
    (int, [float,float], float)
        compid :        The id wrt comps of the next component of polyline B, 
                        or -1 if no next component is found
        [float,float] :  The updated effective interval of component
        float :          The distance from the switch point between next
                         and previous component to segment a
        float :          The k-value of the switch point
    """
    # initialize variables
    new_comp_id = -1
    new_ei = []
    k_min = 2
    # check all components, looking for one with earliest switch point with previous component
    for cand_id in range(len(comps)):
        if cand_id != prev_id: # can't return same component!
            cand_dr = drs[cand_id] #  distance representation of new candidate
            cand_ei = eis[cand_id] #  effective interval of new candidate
            if len(cand_ei) > 0: # check if component has an effective interval
                ks = hu.switchPoint(drs[prev_id], cand_dr) # k-values of switch points between current component and new candidate
                for k in ks:
                    if 0 < k < 1: # switch point within segment
                        if prev_ei[0] <= k <= prev_ei[1]:  # switch point in effective interval of previous component
                            if cand_ei[0] <= k <= cand_ei[1]:   # switch point in effective interval of new candidate
                                if cand_ei[1] > prev_ei[0]: # moving forward    
                                    if k < k_min: # we've got a new winner
                                        k_min = k
                                        new_comp_id = cand_id
                                        new_ei = (k_min,cand_ei[1])
    # If no switch point is found and we haven't reached the end of the segment, 
    # move to next node or vertex
    # *** Hangout says to move to next sequential vertex on B, but I think we have to check both
    #         sides
    if new_comp_id == -1:
        if prev_ei[1] < 1:
            cand_comps = []
            b = comps[prev_id][1]
            if comps[prev_id][0] == True: # segment
                cand_comps = [(False,b+i) for i in [0,1] if hu.compValid(len(B),(False,b+i))]
            else: # vertex
                cand_comps = [(True,b+i) for i in [-1,0] if hu.compValid(len(B),(True,b+i))]
            cand_ids = [comps.index(c) for c in cand_comps if c in comps]
            # choose component with higher effective interval    
            if len(cand_ids) > 0:
                new_comp_id = max(cand_ids, key = lambda id: max(eis[id]))
                new_ei = hu.withinUnitInterval(eis[new_comp_id][0],eis[new_comp_id][1])
                # use new component only if we are moving forward along a
                if max(new_ei) < min(prev_ei):
                    return (-1,None,None,None)
    # check that a new component was found
    if new_comp_id == -1:
        return (-1,None,None,None)
    else:
        # calculate distance at switch point        
        k = new_ei[0]        
        len_a = g.distance(A[a],A[a+1])
        d = hu.componentDistance(drs[new_comp_id], k, len_a)
        return (new_comp_id,new_ei,d,k)
