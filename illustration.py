# -*- coding: utf-8 -*-
"""
***********************************************
                 BASIC USAGE
***********************************************
"""

# Import the main module
import polyline_hausdorff as ph

# Create two polylines
A = [(15, 1), (28, 11), (13, 26), (1, 18)]
B = [(26, 8), (13, -2), (13, 19), (10, 19), (10, 23), (1, 17)]

# Compute the hausdorff distance. This will return a tuple of three items.
hausdorff_info = ph.polyline_hausdorff(A,B)

# The hausdorff distance is the first item in the tuple
h_dist = hausdorff_info[0]
print("The true Polyline Hausdorff distance is {:.3f}.".format(h_dist))




"""
***********************************************
                 PLOT RESULTS
***********************************************
"""

import matplotlib.pyplot as plt
import math

# plot the original polylines
x,y = zip(*A)
plt.plot(x,y,linewidth=3,c='darkorange',marker='o',markersize=5,mfc='white',mec='gray')
x,y = zip(*B)
plt.plot(x,y,linewidth=3,c='blue',marker='o',markersize=5,mfc='white',mec='gray')
plt.axis('equal')

# plot the source location
srcloc = hausdorff_info[1]  
plt.plot(srcloc[0],srcloc[1],marker='X',mec='red',mfc='red',markersize=12)

# plotting text and arrows is complicated - let's make a function for it
def showDist(srcloc,trgloc,lbl,augment=0):
    # draw arrow from source to target
    plt.arrow(srcloc[0],srcloc[1],(trgloc[0]-srcloc[0]),(trgloc[1]-srcloc[1]),color='black',head_width=0.25,head_length=0.25,overhang=0.15,length_includes_head = True)
    # calculate rotation parameters for text
    r = math.degrees(math.atan2(trgloc[1]-srcloc[1],trgloc[0]-srcloc[0]))% 360
    if 100 < r < 260:
        r -= 180
    # calculate text placement
    offset = 0.5
    textx,texty = -1*offset*math.sin(math.radians(r)) + (srcloc[0]+trgloc[0])/2, offset*math.cos(math.radians(r)) + (srcloc[1]+trgloc[1])/2    
    # plot text
    plt.text(textx,texty,lbl,ha='center',va='center',fontsize=15+augment,rotation=r)
    # calculate distance
    dx2 = (trgloc[0] - srcloc[0])**2 
    dy2 = (trgloc[1] - srcloc[1])**2 
    d = math.sqrt(dx2+dy2)
    # show on plot
    offset=-0.5
    textx,texty = -1*offset*math.sin(math.radians(r)) + (srcloc[0]+trgloc[0])/2, offset*math.cos(math.radians(r)) + (srcloc[1]+trgloc[1])/2    
    plt.text(textx,texty,"{:.2f}".format(d),ha='center',va='center',fontsize=12+augment,rotation=r)    

# plot the target locations, with text and arrows
trglocs = hausdorff_info[2]
for trgloc in trglocs:
    showDist(srcloc,trgloc,"polyline hausdorff")




"""
***********************************************
COMPARE WITH SCIPY
***********************************************
"""
try:
    from scipy.spatial.distance import directed_hausdorff
    import numpy as np
    # create numpy arrays
    u = np.array(A)
    v = np.array(B)
    # calculate hausdorff distance and get participating vertices
    d1,a1,b1 = directed_hausdorff(u,v)
    d2,b2,a2 = directed_hausdorff(v,u)
    if d1 > d2:
        d = d1
        srcloc = A[a1]
        trgloc = B[b1]
    else:
        d = d2
        srcloc = B[b2]
        trgloc = A[a2]
    # print results
    print("SciPy Hausdorff distance is {:.2f}.".format(d))
    # show graphically
    showDist(srcloc,trgloc,"scipy")
except ImportError:
    print("Install scipy & numpy to see a comparison with SciPy's Directed Hausdorff distance.")




"""
***********************************************
COMPARE WITH SHAPELY
***********************************************
"""
import utils_hausdorff as uh
# shapely appears to calculate the largest vertex-to-edge gap
try:
    
    from shapely.geometry import LineString
    # create linestring objects
    lineA = LineString(A)
    lineB = LineString(B)
    # calculate hausdorff distance
    shapely_hausdorff = lineA.hausdorff_distance(lineB)
    # print results
    print("Shapely Hausdorff distance is {:.2f}.".format(shapely_hausdorff))
    # show graphically
    srcloc = B[2] 
    trgloc = uh.nearLoc(srcloc, A, (True,1)) 
    showDist(srcloc,trgloc,"shapely")
except ImportError:
    print("Install shapely to see a comparison with shapely's Hausdorff distance.")



    
plt.show()
print("finished.")