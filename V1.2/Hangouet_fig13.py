# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 23:11:10 2021

@author: bjkronenfeld
"""

import matplotlib.pyplot as plt
import math
import numpy as np
import utils_geom as geom
import utils_hausdorff as hu
import matplotlib.patheffects as pe
      
def cbsafe(n):
    """
    Color-blind safe colors (1-9)
    """    
    CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
    return CB_color_cycle[n%9]

# input lines from Hangouet 1995, fig. 13
# slightly modified so that solution is as indicated by Hangouet
A = [(0, 0), (5, 0)]
B = [(2.4, -1), (2.5, -0.4), (4.5, -0.9), (6.4, 0.8), (4, 0.5), (2, 2.3), (0.5, 0.5)]


# set up plot
plt.figure(figsize=(6,4),dpi=128)
# plot the original polylines
x,y = zip(*A)
plt.plot(x,y,linewidth=2,c='darkorange',marker='o',markersize=8,mfc='white',mec='darkorange')
x,y = zip(*B)
plt.plot(x,y,linewidth=2,c='mediumblue',marker='o', markersize=8,mfc='white',mec='mediumblue')
plt.axis('equal')

# label vertices
i=0
v_offset = [-1,1,-1,1,-1,1,-1]
offd = 0.23
for bx, by in B:
    i+=1
    plt.text(bx,by+v_offset[i-1]*offd,"{}".format(i),ha='center',va='center', color="mediumblue", path_effects=[pe.withStroke(linewidth=2, foreground="white")])

# plot distances to vertices
bs = [6,4,1]
for b in bs:
    vdr = hu.vertDistRep(A, B, 0, b)
    vei = hu.vertEffectiveInterval(A, B, 0, b)
    parab = hu.vertDistAcross(vdr,vei,5)
    parab = [(x[0]*5,x[1]) for x in parab]
    parab = [x for x in parab if x[1] < 3]
    x,y = zip(*parab)
    plt.plot(x,y,'--',linewidth=0.8,c='black',marker='o',markersize=0,mfc='white',mec='gray')

# plot distances to segments
for s in [0,1,2,3,4,5]:
    # get distance representation
    sdr = hu.segDistRep(A,B,0,s)
    sei = hu.segEffectiveInterval(A, B, 0, s)
    dist_line = hu.segDistAcross(sdr,sei,5)
    dist_line = [x for x in dist_line if x[1] < 3]
    x,y = zip(*dist_line)
    x = [a*5 for a in x]
    plt.plot(x,y,'--',linewidth=0.8,c='black',marker='o',markersize=0,mfc='white',mec='gray')

# set up axes
plt.xlim(-2,8)    
plt.ylim(-1.5,4)