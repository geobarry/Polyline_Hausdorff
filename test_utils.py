# -*- coding: utf-8 -*-
"""
test data for testing and illustrating Hausdorff distance
"""

import matplotlib.pyplot as plt
import math

def case(casenum):
    """
    Creates a test case consisting of two polylines

    Parameters
    ----------
    casenum : int
        Number between 0 and count -1.

    Returns
    -------
    Two polylines, or the count if casenum == -1

    """
    cases = []
    
    # input lines from Hangouet 1995, fig. 13
    # slightly modified so that solution is as indicated by Hangouet
    A = [(5,0.95),(0.44,-1.99),(0,0)]
    B = [(3.5,-1.2),(3.19,-0.68),(5,0.3),(5.61,1.77),(3.19,0.59),(0.18,0.59),(0.77,-1.6)]
    cases.append((A,B))
    
    # parallel lines
    A = [(0,0),(3,3),(12,2),(14,5)]
    B = [(0,12),(4,6),(13,5),(9,6)]
    cases.append((A,B))
    
    # coincident lines
    A = [(0,0),(4,6),(12,2),(14,5),(17,7)]
    B = [(0,10),(4,6),(12,2),(9,7)]
    cases.append((A,B))
    
    # overlapping lines
    A = [(0,0),(4,6),(12,2),(14,5),(17,7)]
    B = [(0,10),(6,5),(10,3),(9,7)]        
    cases.append((A,B))

    # two parallel segments in a row
    A = [(0,0),(4,6),(14,1),(16,3),(19,5)]
    B = [(0,10),(6,6),(8,5),(10,4),(9,7)]        
    cases.append((A,B))
    
    # two overlapping segments in a row
    A = [(0,0),(4,6),(14,1),(16,3),(19,5)]
    B = [(0,10),(6,5),(8,4),(10,3),(9,7)]        
    cases.append((A,B))
    
    # * cleanup *
    if casenum == -1:
        return len(cases)
    else:
        return cases[casenum]
    
def plot_hausdorff_solution(A,B,d,loc,title,A_highlights = [], B_highlights = []):
    # Plots the input lines with the location (of the far Hausdorff location) 
    # and a circle at distance d from the location
    # If highlight components are provided, these will be highlighted
    # plot input lines
    ax,ay=zip(*A)
    bx,by=zip(*B)
    plt.plot(ax,ay, color = "orange")
    plt.plot(bx,by, color = "blue")
    plt.axis('equal')
    off = (-0.1,-0.1)
    for i in range(len(A)):
        a = A[i]
        plt.text(a[0]+off[0],a[1]+off[1],i)
    off = (0.1,0.1)
    for i in range(len(B)):
        b = B[i]
        plt.text(b[0]+off[0],b[1]+off[1],i)
    # highlighted segment
    for comp in A_highlights:
        if comp[0]: # segment
            s = comp[1]
            highlightx=[A[s][0],A[s+1][0]]
            highlighty=[A[s][1],A[s+1][1]]
            plt.plot(highlightx,highlighty,color='darkorange',linewidth=3)
        else: # vertex
            v = comp[1]
            plt.plot(A[v][0],A[v][1],marker = 'o', mec = 'darkorange')
    for comp in B_highlights:
        if comp[0]: # segment
            s = comp[1]
            highlightx=[B[s][0],B[s+1][0]]
            highlighty=[B[s][1],B[s+1][1]]
            plt.plot(highlightx,highlighty,color='darkblue',linewidth=3)
        else: # vertex
            v = comp[1]
            plt.plot(B[v][0],B[v][1],marker = 'o', mec = 'darkblue')
    # solution point
    plt.plot(loc[0],loc[1],'ro')
    # circle
    degs = list(range(361))
    rads = [math.radians(x) for x in degs]
    pts = [(loc[0]+d*math.sin(r),loc[1]+d*math.cos(r)) for r in rads]
    x,y = zip(*pts)
    plt.plot(x,y,color='grey')
    pretitle = "A: orAnge     B: Blue\n"
    plt.suptitle(pretitle + title)
    plt.show()