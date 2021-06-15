# -*- coding: utf-8 -*-
"""
Created on Wed May 26 13:29:02 2021

@author: bjkronenfeld
"""
import math as m
import utils_geom as g

def componentDistance (dr,k,len_a):
    """
    Calculates the distance from a location on segment a on polyline A to 
    a component of polyline B.
    
    Parameters
    ----------
    A : list of (x,y) tuples
        First polyline.
    dr : distance representation
        distance representation of component of polyline B.
    k : float
        k-value of location along segment a.
    a : int
        ID of segment on polyline A.

    Returns
    -------
    Float
        Distance from location k on segment a to dr component.
    
    # LJ and MT
    """
    
    if dr[0] == True:
        answer = segDistance(dr,k,len_a)
    else:    
        answer = vertDistance(dr,k,len_a)
    return answer

def component_label(comp):
    """for debugging"""
    if comp == None:
        return "None"
    else:
        if comp[0] == True:
            return"seg {}".format(comp[1])
        else:
            return "vert {}".format(comp[1])

def compValid(n_vert,comp):
    """
    Determines if the input component is valid on the input polyline.
    Args:
        polyline : list of tuple
            The input polyline.
        comp : (bool, float)
            The input component
        
    Returns:
        Bool
    
    # BK
    """
    if comp[1] < 0:
        return False
    if comp[0]: # component is a segment
        if comp[1] > n_vert-2:
            return False
    else: # component is a vertex
        if comp[1] > n_vert-1:
            return False
    return True

def distanceRepresentation (A,B,a,bcomp):
    """
    Computes the distance representation of a component of B with respect to segment a 
    Args:
        A,B - the two polylines
        a (int) - segment of A
        bcomp  - componet of B 
        
    Returns:
        distance representation
    
    # LJ and TJ
    """
    # if bcomp is a vertex
    if bcomp[0] == False:
        answer = vertDistRep(A, B, a, bcomp[1])
        return answer
    #if bcomp is a segment
    else:
        answer = segDistRep(A, B, a, bcomp[1])
        return answer


def effectiveInterval(A,B,a,bcomp):
    """
    Determines the range along segment a that is closer to the component bcomp 
    than to either adjacent component.
    Parameters
    ----------
    A : [(x,y),(x,y)...]
        The main polyline.
    B : [(x,y),(x,y)...]
        The other polyline
    a : Integer
        Segment on A
    bcomp : Component
        component of B

    Returns
    -------
    tuple of two k-values

    @author: Folaitan & TJones
    """
    # if statement to determine if bcomp is a segment or vertex
    if bcomp[0] == True:
        sint = segEffectiveInterval(A, B, a, bcomp[1])
        return sint
    elif bcomp[0] == False:
        vint = vertEffectiveInterval(A, B, a, bcomp[1])
        return vint


def segDistance(dr,k,len_a):
    """
    Computes the distance from dr to k on a segment a.
    
    Parameters
    ----------
    dr : (is_seg,sk or None,sin_theta or q-value)
        dist rep of segament B.
        
    k: k-value along segment a.
    
    len_a: length of segment a. 
    Returns
    ----------
    float
        The distance from dr to k on segment a.
        
    @author: megshithakur and Farouk
    """
    # check if segments are parallel
    if dr[1] == None:
        return dr[2] * len_a
    else: # normal case
        # segDistance formula
        segDist=abs((k-dr[1])*len_a*dr[2])
        return segDist

def segDistRep (A,B,a,b):
    """
    Constructs the component distance representation for segment b with respect to segment a. 
    
    Parameters
    ----------
        A : list of (x,y) tuples
            The first polyline.
        B : list of (x,y) tuples
            The second polyline.
        a : int
            Index of segment on A
        b : int
            Index of segment on B
    
    Returns
    ----------
    tuple (distance representation):
        isSeg : boolean
            True (always True by definition)
        k : float
            k-value of intersection between two segments, or None if segments
            are parallel
        sin_theta : float
            sine of angle between lines through segments, or q-distance between
            lines if they are parallel.
    
    # LJ and TJ
    """
    # initiate list for results
    dist_rep = []  
    # Create a boolean statement saying that this is indeed a segment
    dist_rep.append(True)
    # run intersection tool to find the point of intersection between two infinite lines
    x = g.intersection(A[a], A[a+1], B[b], B[b+1])
    # if we've found an intersection, determine sin theta
    if x[0] != None:
        # If there is an intersection, where is it on the x-axis, and then append it
        K = g.kvalue(x, A[a], A[a+1]) # convert to a k-value
        dist_rep.append(K)
        # Find the angle created by the two lines
        s_rad = g.angle([A[a],A[a+1]],[B[b],B[b+1]])
        # Calculate the sine of the angle
        s = m.sin(s_rad)
        #Append the sine of the angle
        dist_rep.append(s)
    # Handle case of no intersection
    # If sin theta is zero, there really isn't an intersection (but this 
    # might not be caught above due to floating point precision errors)
    if x[0] == None or dist_rep[2] == 0:
        # reset
        dist_rep=[]
        dist_rep.append(True)
        # Append a no result (no intersection) into the final list
        dist_rep.append(None)
        # Find the distance between the two lines
        q = g.distance_to_line(B[b],A[a],A[a+1])
        # normalize by length of a
        q=q/g.distance(A[a], A[a+1])
        # Append the distance
        dist_rep.append(q)
    dist_rept = tuple(dist_rep)
    return dist_rept


def withinUnitInterval(a,b):
    """
    Computes the portion of the input interval that is within the interval [0,1]

    Parameters
    ----------
    interval : (float,float)
        An effective interval, not necessarily in sequence.

    Returns
    -------
    The portion of the input interval within [0,1], in sequence from low to high,
    or (-inf,-inf) if the input interval does not overlap the unit interval

    """
    kmin = min(min(a,b),1)
    kmax = max(max(a,b),0)
    if kmin==1 or kmax==0:
        return (float('-inf'),float('-inf'))
    
    else:
        return (kmin,kmax)

def segEffectiveInterval(A,B,a,b, tolerance=0.000001):
    """
    Computes the effective interval of segment b on segment a, that is the interval
    on a for which the interior of segment b is closer than either endpoint
    
    Parameters
    ----------
        A : list of (x,y) tuples
            The first polyline.
        B : list of (x,y) tuples
            The second polyline.
        a : int
            Index of segment on A.
        b : int
            Index of segment on B.
        tolerance : float
            If seg b's k-values on seg a are within this tolerance, 
            the segments will be treated as perpendicular
    Returns:
        tuple
        k1 : float
            First k-value of effective interval. 
        k2 : float
            Second k-value of effective interval.
            
    @author: megshithakur and tannerjones
    """
    a1 = A[a]
    a2 = A[a+1]
    b1 = B[b]
    b2 = B[b+1]
    # check for perpendicular segments
    prjout1 = g.project_pt_to_line(b1, a1, a2)
    prjout2 = g.project_pt_to_line(b2, a1, a2)
    k1 = g.kvalue(prjout1,a1,a2)
    k2 = g.kvalue(prjout2,a1,a2)
    if abs(k1-k2) <= tolerance:
        # segment is perpendicular, so effective interval is entire segment
        return [0,1]
    else:
        # project out from each b vertex to segment a 
        k1 = g.project_out(a1,a2,b1,b2)
        k2 = g.project_out(a1,a2,b2,b1)
        # return them in sequence, bound to range [0,1]
        return withinUnitInterval(k1,k2)

def switchPoint (dr1, dr2):
    """
    Determines point along a where the nearest point on B switches from the 
    component represented by dr1 to the component represented by dr2.

    Parameters
    ----------
    dr1 : dr1
        Distance representation 1
    dr2 : dr2
        Distance representation 2

    Returns
    -------
    [float]	list of k-values

    """
    # input: segment, segment - return Seg Seg switch point
    if dr1[0] == True and dr2[0] == True:
        segseg = segSegSwitchPoint(dr1, dr2)
        return segseg
    # input: vertex, segment - return vert seg switch point
    elif dr1[0] == False and dr2[0] == True:
        vertseg = vertSegSwitchPoint(dr1, dr2)
        return vertseg
    # input: segment, vertex - return vert seg switch point
    elif dr1[0] == True and dr2[0] == False:
        vertseg = vertSegSwitchPoint(dr2, dr1)
        return vertseg
    # input: vertex, vertex - return vert vert switch point
    elif dr1[0] == False and dr2[0] == False:
        vertvert = vertVertSwitchPoint(dr1, dr2)
        return vertvert


def vertDistance(dr, k, len_a):
    """
    Computes the distance from the vertex represented by dr to the location k 
    on the segment a that the distance representation was constructed from.
    Parameters
    ----------
    dr : (is_seg, vk,q)
        distance representation of a vertex of B, being a tuple of three values
    k : the k-value along segment a
    len_a : the length of segment a

    Returns
    -------
    float
        the distance from point k on segment a to the vertex of B represented 
        by dr
        
    @author: Folaitan & @Ljansen
    """
    d = ((k-dr[1])**2) + (dr[2]**2)
    d = m.sqrt(d)
    d = len_a * d
    return d

def vertDistRep(A,B,a,b):
    """
    Constructs the Distance Representation for vertex b with respect to segment a.
    
    Parameters
    ----------
        A : list of (x,y) tuples
            The first polyline.
        B : list of (x,y) tuples
            The second polyline.
        a : int
            Index of segment on A
        b : int
            Index of vertex on B
    Returns
    ----------
    tuple (distance representation): 
        is_seg : bool 
            False(always False by definition).
        k : float 
            k-value of the location of the perpendicular projection of b onto the line through a.
        q : float 
            q-distance from b to the line through a.
            
    @author: MT and TJ
    """
    # initiate list for results
    FinalList = []
    #Create a boolean statement saying that this is not a segment
    FinalList.append(False)
    #run project to line tool to find perpendicular intersection point between point b and segment a 
    p = g.project_pt_to_line(B[b], A[a], A[a+1])
    k = g.kvalue(p, A[a], A[a+1])
    #append result to results list
    FinalList.append(k)
    #run distance tool to find distance between the new point (k) and b
    q = g.distance(p, B[b])
    q = q/g.distance(A[a],A[a+1])
    #append the result
    FinalList.append(q)
    FinalList = tuple(FinalList) 
    return FinalList

def vertEffectiveInterval(A,B,a,b,tolerance=0.000001):
    """
    Computes the effective interval of vertex b on segment a, that is the interval
    on a for which vertex b is closer than either adjacent segment
    
    Parameters
    ----------
        A : list of (x,y) tuples
            The first polyline.
        B : list of (x,y) tuples
            The second polyline.
        a : int
            Index of segment on A
        b : int
            Index of vertex on B
        tolerance : float
            If b segs' k-values on seg a are within this tolerance, 
            the segments will be treated as perpendicular
    Returns:
    ----------
        tuple of two floats, with (-inf,-inf) representing no effective
        interval.
        
    
    @author: megshithakur and tannerJones
    """
    # get coordinates
    a1, a2 = A[a],A[a+1]

    
    # handle cases where b is an end of the polyline
    if b==0 or b==len(B)-1:
        # get neighbor segment id
        if b==0:
            nb = 1
        else:
            nb = len(B)-2
        # get coordinates of b vertex and neighbor
        bc,bnb = B[b],B[nb] 
        
        # determine projections onto segment a
        prjc = g.project_pt_to_line(bc,a1,a2)
        prjnb = g.project_pt_to_line(bnb,a1,a2)
        # determine k-values
        kc = g.kvalue(prjc, a1, a2)
        knb = g.kvalue(prjnb, a1, a2)

        # check for perpendicularity
        if abs(kc-knb) <= tolerance:
            # calculate distances from each B vertex to segment a
            distc = g.distance(bc,prjc)
            distnb = g.distance(bnb,prjnb)
            # is vertex b on segment a?
            if distc < tolerance:
                return (0,1)
            # is the neighboring vertex on segment a?
            elif distnb < tolerance:
                return (float('-inf'),float('-inf'))
            # are the two vertices on opposite sides of segment a?
            else: 
                # calculate areas
                areac = g.area([a1,a2,bc,a1])
                areanb = g.area([a1,a2,bnb,a1])
                # are the two vertices on opposite sides of segment a?
                if (areac > 0) != (areanb > 0):
                    return (float('-inf'),float('-inf'))
                # is the center vertex closer to segment b than the neighbor
                elif distc < distnb:
                    return [0,1]
                else:
                    return (float('-inf'),float('-inf'))
        else:
            # project points out from b to segment a
            kbout = g.project_out(a1, a2, bc, bnb)
            # return interval from kbout to end of segment a opposite neighbor
            if knb < kc:
                return withinUnitInterval(kbout,1)
            else:
                return withinUnitInterval(0,kbout)
    else: # vertex b has two neighbors     
        # get coordinates of B vertices
        bp,bc,bn = B[b-1],B[b],B[b+1] # prev, current, next vertices on B
        # project each vertex of B onto segment A
        prjp = g.project_pt_to_line(bp,a1,a2)
        prjc = g.project_pt_to_line(bc,a1,a2)
        prjn = g.project_pt_to_line(bn,a1,a2)
        # get k-values of projections of B vertices onto a
        kp = g.kvalue(prjp, a1, a2)
        kc = g.kvalue(prjc, a1, a2)
        kn = g.kvalue(prjn, a1, a2)    
        # check for perpendicular segments
        prev_perp = abs(kc-kp) <= tolerance
        next_perp = abs(kn-kp) <= tolerance
        if prev_perp and next_perp: # both perpendicular
            # get distances and areas to all vertices
            dcur = g.distance(bc,prjc)
            dprev = g.distance(bp,prjp)
            dnext = g.distance(bp,prjn)
            areacur = g.area([a1,a2,bc,a1])
            areaprev = g.area([a1,a2,bp,a1])
            areanext = g.area([a1,a2,bn,a1])
            # vertex has effective interval only if it is closest to segment a
            # and all vertices are on same side of segment a
            if (areaprev > 0) == (areanext > 0) and (areacur > 0) == (areanext > 0) and dcur < dprev and dcur < dnext:
                return (0,-1)
            else:
                return (float('-inf'),float('-inf'))
        elif prev_perp or next_perp: # one perpendicular
            # get coordinates and k-values of vertex on perpendicular segment, other vertex
            if prev_perp:
                bperp = bp
                bother = bn
                kother = kn
                prjperp = prjp
            else:
                bperp = bn
                bother = bp
                kother = kp
                prjperp = prjn            
            # get distances and areas of each vertex on perpendicular segment
            dperp = g.distance(bperp,prjperp)
            dcur = g.distance(bc,prjc)
            areaperp = g.area([a1,a2,bperp,a1])                    
            areacur = g.area([a1,a2,bc,a1])                    
            # Is perpendicular segment on one side of A and current vertex is closer?
            if (areaperp > 0) == (areacur > 0) and dcur < dperp:
                # project out from other segment
                otherprjoutk = g.project_out(a1,a2,bc,bother)
                # interval is from end of line opposite other to k-value of other
                if kother > kc:
                    return (0,otherprjoutk)
                else:
                    return (otherprjoutk,1)
            else:
                return (float('-inf'),float('-inf'))
        else: # neither perpendicular
            # get k-values of projections of vertex out from each B segment onto segment a            
            kcpout = g.project_out(a1, a2, bc, bp) 
            kcnout = g.project_out(a1, a2, bc, bn)
            # check sides of neighboring vertices with respect to b
            if kp < kc and kn < kc: # both neighbors left of b
                maxk = max(kcpout,kcnout)
                return withinUnitInterval(maxk,1)
            elif kp > kc and kn > kc: # both neighbors right of b            
                mink = min(kcpout,kcnout)
                return withinUnitInterval(0,mink)
            else: # nieghbors on either side of b
                # determine min and max k-values of b based on positions of
                # previous and next vertices
                if kp < kc: # previous neighbor left
                    mink = kcpout
                    maxk = kcnout
                else: # previous neighbor right
                    mink = kcnout
                    maxk = kcpout
                if mink <= maxk: 
                    return withinUnitInterval(mink,maxk)
                else:
                    return (float('-inf'),float('-inf'))

def segSegSwitchPoint(dr1, dr2):
    """
    Determines the k-values of the two locations along segment a that are equidistant 
    to the two segments b1 and b2 represented by dr1 and dr2, i.e. the location at 
    which the nearest component of a “switches” from b1 to b2.

    Parameters
    ----------
    dr1 : seg_p1
        distance representation of a segment on b 
    dr2 : seg_p2
        distance representation of another segment on b 

    Returns
    -------
    List of 1 or 2 floats representing k-value(s) of the switch points
    
    @author: Folaitan & Ljansen
    """
    k = []
    if dr1[1] == None and dr2[1] == None: # both segments parallel to a
        pass
    elif dr1[1] == None: # first segment parallel to a
        k_out_1 = dr2[1] + (dr1[2]/dr2[2])
        k_out_2 = dr2[1] - (dr1[2]/dr2[2])
        k.append(k_out_1)
        k.append(k_out_2)
    elif dr2[1] == None: # second segment parallel to a
        k_out_1 = dr1[1] + ((dr2[2])/dr1[2])
        k_out_2 = dr1[1] - ((dr2[2])/dr1[2])
        k.append(k_out_1)
        k.append(k_out_2)
    else: # neither segment parallel to a
        # create more readable variables
        k1 = dr1[1]
        k2 = dr2[1]
        # calculate alpha parameter
        a = dr2[2]/dr1[2]
        if a == 1: # first solution is not valid; return second solution
            k_out_2 = (k1+a*k2)/(1+a)
            k.append(k_out_2)
        elif a == -1: # second solution is not valid; return first solution
            k_out_1 = (k1-a*k2)/(1-a)
            k.append(k_out_1)
        else: # return both solutions
            k_out_1 = (k1-a*k2)/(1-a)
            k.append(k_out_1)
            k_out_2 = (k1+a*k2)/(1+a)
            k.append(k_out_2)
    # return values in ascending order, for consistency
    return sorted(k)

def vertSegSwitchPoint(vdr, sdr):
    """
    Determines the k-value of the location along segment a that is equidistant 
    to the input vertex and input segment, i.e. the location at which the nearest 
    component “switches” from the first component to the second.
    Parameters
    ----------
    vdr : ver_rep : (False, k, q)
        distance representation of a vertex on B
    sdr : seg_rep : (True, k, sin_theta) or (True, none, q)
        distance representation of a segment on B

    Returns
    -------
    [float,float]
    List containing two floats representing the k-values of the switch points
    
    @author: Folaitan & Ljansen
    """
    r_list = []
    # get distance representation values into more readable variables
    k_vert = vdr[1]
    q_vert = vdr[2]
    k_seg = sdr[1]
    sin_theta = sdr[2]
    # three cases
    if sin_theta== 1: # b segment is perpendicular to A
        numerator = ((k_seg**2)-(k_vert**2)-(q_vert**2))
        denominator = (2*k_seg) - (2*k_vert)
        if denominator == 0: 
            # point and line have same k-value, so this is the switch point
            return [k_seg]
        vertSegPoint = numerator/denominator
        r_list.append(vertSegPoint)
    else:
        if sdr[1] == None: # b segment is parallel to a
            q_seg = sin_theta # distance representation value is q not sin_theta
            a = 1
            b = -2*k_vert
            c = k_vert**2 + q_vert**2 - q_seg**2
        else: # normal case
            a = (sin_theta**2)-1
            b = (2*k_vert)-((2*k_seg)*(sin_theta**2))
            c = ((k_seg**2)*(sin_theta**2))-(k_vert**2)-(q_vert**2)
        inside_root = (b**2)-(4*a*c)
        # catch floating point precision issues:
        # better to show a switch point when there is none than to miss one
        if -0.0000000000001 < inside_root < 0:
            inside_root = 0                
        # if inside_root is less than zero, quadratic formula has no 
        # solution and there is no switch point to return
        # if it equals zero, there is one switch point
        if inside_root >= 0:
            r_list.append((-b + (m.sqrt(inside_root)))/(2*a))
        if inside_root > 0:
            r_list.append((-b - (m.sqrt(inside_root)))/(2*a))                
    # let's put these in ascending order, just to be safe
    return sorted(r_list)

            
def vertVertSwitchPoint (dr1,dr2):
    """
    Determines the k-value of the location along segment a that is equidistant 
    to the two vertices b1 and b2 represented by dr1 and dr2, i.e. the location 
    at which the nearest component of a “switches” from b1 to b2. 

    Parameters
    ----------
    dr1 : tuple, distance representation of a vertex of B
    dr2 : tuple, distance representation of a different vertex of B

    Returns : [k] - list containing one float representing the k-value of the switch point
    
    # MT and LJ
    """
    #checks to see if k values for dr1 and dr2 are not equal
    if dr1[1] != dr2[1]:
        #calculates k value for switch point
        k = [((dr2[2]**2 - dr1[2]**2) + (dr2[1]**2 - dr1[1]**2)) / (2 * (dr2[1]-dr1[1]))] 
    else:
        #k is an empty list
        k = []
    return k

def candidateComponents(A,B,a):
    """
    Identifies all components of B that could be the target of the Hausdorff 
    distance from segment a on A.

    Parameters
    ----------
    A : [(x,y),...] list of tuples 
        The coordinates of the main polyline.
    B : [(x,y),...] list of tuples 
        The coordinates of the other polyline.
    a : int
        The index of a segment on polyline A.
    Returns
    ----------
    [(bool, int)]
        A list of candidate components on B.
    """
    # for now, simply return a list of all components of B
    result = []
    # get vertex components
    for i in range(len(B)):
        result.append((False,i))
    for i in range(len(B)-1):
        result.append((True,i))
    return result

def nearSegment(A,B,a):
    """
    Identifies the segment of of B nearest to vertex a on A.

    Parameters
    ----------
    A : [(x,y),...] list of tuples 
        The coordinates of the main polyline.
    B : [(x,y),...] list of tuples 
        The coordinates of the other polyline.
    a : int
        The index of a vertex on polyline A.
    Returns
    ----------
    int
        The index of the nearest segment on B.
    """
    # for now, this will be coded with a "brute force" method, checking the 
    # distance from every segment of B to vertex a of A
    # later this should be updated to use an indexing structure such as 
    # an r-tree for computational efficiency
    # initialize to first segment
    min_index = 0
    min_d = g.distance_to_segment(A[a], B[0], B[1])
    # check other segments
    for b in range(1,len(B)-1):
        d = g.distance_to_segment(A[a],B[b],B[b+1])
        if d < min_d:
            min_d = d
            min_index = b
    return min_index

def nearComponent(A,B,a,b):
    """
    Among segment b, vertex b and vertex b+1, determines which component is nearest to vertex a.

    Parameters
    ----------
    A : [(x,y),...] list of tuples 
        The coordinates of the main polyline.
    B : [(x,y),...] list of tuples 
        The coordinates of the other polyline.
    a : int
        The index of a vertex on polyline A.
    b : int
        The index of a segment on polyline B.
    Returns
    ----------
    (component, float)
        The nearest component of B along with its distance from vertex a.
    """
    # # get distances from each component
    # d_seg = g.distance_to_segment(A[a],B[b],B[b+1])
    # d_vert1 = g.distance(A[a],B[b])
    # d_vert2 = g.distance(A[a], B[b+1])
    # # if distance to d_seg is strictly the minimum, return the segment
    # if d_seg == min(d_seg,d_vert1,d_vert2):
    #     return ((True,b),d_seg)
    # elif d_vert1 < d_vert2:
    #     return((False,b),d_vert1)
    # else:
    #     return ((False,b+1),d_vert2)
    
    # calculate k-value of projection of vertex a onto segment b
    prj = g.project_pt_to_line(A[a], B[b], B[b+1])
    k = g.kvalue(prj, B[b], B[b+1])
    if k <= 0: # return vertex b
        d = g.distance(A[a],B[b])
        return ((False,b),d)
    elif k >= 1: # return vertex b+1
        d = g.distance(A[a],B[b+1])
        return ((False,b+1),d)
    else: # return segment b
        d = g.distance(A[a],prj)
        return ((True,b),d)        

def nearLoc(srcloc,trgline,trgcomp):
    """
    Finds the location on the target component nearest to the source
    location.

    Parameters
    ----------
    srcloc : (float,float)
        The coordinates of the source location.
    trgline : [(x,y),...] list of tuples 
        The coordinates of the target polyline.
    trgcomp : (bool,int)
        The target component.
    
    Returns
    ----------
    (float,float)
        The nearest location on the target component to the source location.
    """
    if trgcomp[0] == False: # target is a vertex
        return trgline[trgcomp[0]]
    else: # target is a segment
        trgstart = trgline[trgcomp[1]]
        trgend = trgline[trgcomp[1]+1]
        srcprj = g.project_pt_to_line(srcloc, trgstart,trgend)
        k = g.kvalue(srcprj,trgstart,trgend)
        if k <= 0:
            return trgstart
        elif k >=1:
            return trgend
        else:
            x = trgstart[0] + k * (trgend[0]-trgstart[0])
            y = trgstart[1] + k * (trgend[1]-trgstart[1])
            return (x,y)
        
def checkSegment(c1,c2):
    """
    Determines whether or not it is necessary to further process a segment
    given that its endpoints have been processed already.

    Parameters
    ----------
    c1 : component
        The component of B closest to the first vertex of a segment of A.
    c2 : component
        The component of B closest to the second vertex of a segment of A.
    Returns
    ----------
    (component, float)
        The nearest component of B along with its distance from vertex a.
    """
    # No need to check a segment if either:
    # c1 and c2 are the same component, or
    # c1 and c2 are consective vertices
    if c1==c2:
        return False
    elif c1[0] == False and c2[0] == False and abs(c1[1]-c2[1])==1:
        return False
    else:
        return True
