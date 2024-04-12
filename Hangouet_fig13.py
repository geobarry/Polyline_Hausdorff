# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 23:11:10 2021

@author: bjkronenfeld
"""

import matplotlib.pyplot as plt
import math
import numpy as np
import polyline_hausdorff as h
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
A = [(0, 0), (5, 0)] # plotting code assumes a horizontal line beginning at (0,0)
B = [(2.4, -1), (2.5, -0.4), (4.5, -0.9), (6.4, 0.8), (4, 0.5), (2, 2.3), (0.5, 0.5)]

# colors
A_color = cbsafe(1)  # "darkorange"
B_color = cbsafe(0)# "tab:blue"

# plotting options
plot_distance_functions = True
bf = {
      "v7":0.15,
      "s6":0.07,
      "v2":0.05,
      "s2":0.03,
      "v5":0.06,
      "s4":0.05,
      "s3":0.035
      } # buffer on either side of distance function section

# folder/file to save figure to (None for inline printing only)
folder = r"C:\CaGIS Board Dropbox\cantaloupe bob\Barry\Research\Projects\polyline difference metrics\Hausdorff\images"
file = "Hangouët_with_distance_functions_transparent.png"

# precalculate some things
xmin = min(x[0] for x in A+B) 
ymin = min(x[1] for x in A+B)
xmax = max(x[0] for x in A+B)
ymax = max(x[1] for x in A+B) 

trg_area = 15
bf_pct = 10
dx = (xmax-xmin) * (100+2*bf_pct)/100
dy = (ymax-ymin) * (100+2*bf_pct)/100
scale_factor = (trg_area / (dx * dy))**0.5
fig_w = dx * scale_factor
fig_h = dy * scale_factor


full_path = f"{folder}\\{file}"

# get sections of polyline A
a_k_d_comp_rep = h.near_components(A, B)

# set up plot
fig, ax = plt.subplots(figsize=(fig_w,fig_h),dpi=95) #128
# plot the original polylines
x,y = zip(*A)
plt.plot(x,y,linewidth=2,c=A_color,marker='o',markersize=8,mfc='white',mec=A_color)
x,y = zip(*B)

plt.plot(x,y,linewidth=2,c=B_color,marker='o', markersize=8,mfc=(0.99,0.99,0.99),mec=B_color)
plt.axis('equal')

# plot distance functions
len_A = (((A[1][0]-A[0][0])**2)+((A[1][1]-A[0][1])**2))**0.5
print(f"length A: {len_A}")
for i in range(len(a_k_d_comp_rep)-1):
    a,k,d,comp,rep = a_k_d_comp_rep[i]
    k_next = a_k_d_comp_rep[i+1][1]
    # buffer section bounds to display distance functions better
    lbl = hu.component_label((comp[0],comp[1]+1))
    k1 = k - bf[lbl]
    k2 = k_next + bf[lbl]
    # trim segments by effective interval
    if comp[0]:
        ei = hu.effectiveInterval(A, B, a, comp)
        if k1 < ei[0]:
            k1 = ei[0]
        if k2 > ei[1]:
            k2 = ei[1]
    # get k-values, distances of distance function
    dist_line = hu.distAcross(comp,rep,(k1,k2),len_A)
    # plot distance functions 
    # assumes that A is horizontal at y=0
    if plot_distance_functions:
        x,y = zip(*dist_line)
        x = [a*len_A for a in x]
        plt.plot(x,y,'--',linewidth=0.8,c='black',marker='o',markersize=0,mfc='white',mec='gray')
        # plot text label for component
        plt.text(len_A*(k+k_next)/2,0.1,hu.component_label((comp[0],comp[1]+1)),
                 ha = "center",
                 fontsize = 8 
                 )
# plot section dividers
if plot_distance_functions:
    for i in range(len(a_k_d_comp_rep)):
        a,k,d,comp,rep = a_k_d_comp_rep[i]
        plt.plot([k*len_A,k*len_A],[0,d],
                     linewidth = 0.8,
                     c = (0.5,0.5,0.5)
                 )
    
# plot vertices on polyline B
x,y = zip(*B)
plt.plot(x,y,linewidth=0,c=B_color,marker='o', markersize=9,mfc='white',mec=B_color)

# label vertices
i=0
x_offset = [-1,-1.3,0,1,0,0,-0.5]
y_offset = [-1,0.2,-1.3,1,-1,1,1]
offd = 0 # 0.23
for bx, by in B:  
    i+=1
    plt.text(bx+0.01,by-0.01,f"{i}",
             ha='center',
             va='center', 
             color=B_color, 
             fontsize = 7,
             path_effects=[pe.withStroke(linewidth=2, foreground="white")])


# set up axes
#plt.xlim(-2,8)    
#plt.ylim(-1.5,4)
ax.set_axis_off()

# save figure
plt.savefig(full_path, dpi = 600, transparent = True)