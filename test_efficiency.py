# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 19:54:22 2024

@author: Administrator
"""
from polyline_hausdorff import directional_polyline_avg_dist as avg_dist
from linesimplify.missouri_river import missouri_river as missouri_river
from linesimplify import apsc
from linesimplify import douglas_peucker as dp
import shapefile
import time
import win32clipboard
import matplotlib.pyplot as plt


def Missouri():
    # If we want to use a small feature
    return missouri_river()

def cannonball():
    # Ok let's look at a bigger feature
    filename = r"C:\CaGIS Board Dropbox\cantaloupe bob\Research\Projects\s3_metric\data\testStreams50_watt.shp"
    sf = shapefile.Reader(filename)
    return sf.shapes()[2].points

def copy2clip(txt):
    # set clipboard data
    win32clipboard.OpenClipboard()
    win32clipboard.EmptyClipboard()
    win32clipboard.SetClipboardText(txt)
    win32clipboard.CloseClipboard()


def time_test_equal_vertex_counts(feat_a_func,feat_b_func,vertex_counts):
    msg = f"n\telapsed\tavg_dist"
    print(msg)
    data = msg
    for n in vertex_counts:    
        # Create two versions of a simplified line
        feat_a = feat_a_func(n)
        feat_b = feat_b_func(n)
        # calculate average polyline distance
        start = time.time() 
        avgD = avg_dist(feat_a,feat_b,brute_force)
        finish = time.time()
        elapsed = finish-start
        msg = f"{n} \t {elapsed} \t {avgD}"
        data = data + "\n" + msg
        print(msg)
    copy2clip(data)
    print("Finished")

def time_test_unequal_vertex_counts(
        feat_a_func,
        feat_b_func,
        vertex_counts,
        brute_force = False):
    msg = "n-m\t" + "\t".join([str(v) for v in vertex_counts])
    print(msg)
    for n in vertex_counts:    
        times = []
        print(f"n: {n}")
        for m in vertex_counts:
            print(f"...m: {m}")
            # Create two versions of a simplified line
            feat_a = feat_a_func(n)
            feat_b = feat_b_func(m)
            # calculate average polyline distance
            start = time.time() 
            avgD = avg_dist(feat_a,feat_b,brute_force)
            finish = time.time()
            elapsed = finish-start
            times.append(elapsed)
        msg += f"\n{n}\t" + "\t".join([str(t) for t in times])
    copy2clip(msg)
    print("Finished")


def zigzag(n,vertical=False):
    pts = [(i/(n-1),i%2) for i in range(n)]
    if vertical:
        pts = [(y,x) for x,y in pts]
    return pts

def test_same_polyline(feat):
    """
    See if it's really giving a nonzero distance for the same polylines compared with each other'
    """
    simp_table=apsc.simplificationTable(feat) # create the simplification table 
    def get_first_feature(n):
        return apsc.simplifiedLine(simp_table,min_pts=n) 
    time_test_equal_vertex_counts(get_first_feature,get_first_feature,vertex_counts)

def time_test_equal(feat,vertex_counts):
    """ Tests how long it takes to compute average distance between 
        two polylines with equal numbers of vertices """
    # precalculate
    simp_table=apsc.simplificationTable(feat) # create the simplification table 
    errors, sorted_errors = dp.get_errors_sortedErrors(feat) # get displacement distance for each point
    def get_first_feature(n):
        return apsc.simplifiedLine(simp_table,min_pts=n) 
    def get_second_feature(n):
        return dp.simplify_by_numPts(feat, n, errors, sorted_errors) # simplify to 10 pts,[]
    time_test_equal_vertex_counts(get_first_feature,get_second_feature,vertex_counts)

def time_test_unequal(feat,vertex_counts,brute_force = False):
    """ Tests how long it takes to compute average distance between 
        two polylines with equal numbers of vertices """
    # precalculate
    simp_table=apsc.simplificationTable(feat) # create the simplification table 
    errors, sorted_errors = dp.get_errors_sortedErrors(feat) # get displacement distance for each point
    def get_first_feature(n):
        return apsc.simplifiedLine(simp_table,min_pts=n) 
    def get_second_feature(n):
        return dp.simplify_by_numPts(feat, n, errors, sorted_errors) # simplify to 10 pts,[]
    time_test_unequal_vertex_counts(
        get_first_feature,
        get_second_feature,
        vertex_counts,
        brute_force = brute_force)


def time_test_zigzag(vertex_counts):
    def zigzag_vertical(n):
        return zigzag(n,False)
    def zigzag_horizontal(n):
        return zigzag(n,True)
    time_test_equal_vertex_counts(zigzag_vertical,zigzag_horizontal,vertex_counts)

def create_example_plots(orig_feat):
    """
    create sample plots for paper demonstrating scenarios with
    differences in computational complexity
    """
    c1 = "blue"
    c2 = "darkorange"
    fig, axes = plt.subplots(1,3,figsize=(7,2))
    plt.axis('equal')
    for ax in axes:
        ax.axis("off")
        plt.axis('equal')
    # show two polylines with equal number of vertices
    n = 25
    simp_table=apsc.simplificationTable(orig_feat) # create the simplification table 
    errors, sorted_errors = dp.get_errors_sortedErrors(orig_feat) # get displacement distance for each point
    feat_a = apsc.simplifiedLine(simp_table,min_pts=n) 
    feat_b = dp.simplify_by_numPts(orig_feat, n, errors, sorted_errors)
    ax = axes[0]
    x,y = zip(*feat_a)
    ax.plot(x,y,color=c1,linewidth=0.5)
    x,y = zip(*feat_b)
    ax.plot(x,y,color=c2,linewidth=0.5)
    ax.set_title('best case:\n equal # vertices', y = -0.3)
    # show one polyline with 10x vertices as the other
    n = 10
    feat_a = apsc.simplifiedLine(simp_table,min_pts=n*10) 
    feat_b = dp.simplify_by_numPts(orig_feat, n, errors, sorted_errors)
    ax = axes[1]
    x,y = zip(*feat_a)
    ax.plot(x,y,color=c1,linewidth=0.5)
    x,y = zip(*feat_b)
    ax.plot(x,y,color=c2,linewidth=0.5)
    ax.set_title('typical case:\n unequal # vertices', y = -0.3)
    # show worst case scenario - 2 zigzags in opposite directions    
    ax = axes[2]
    feat_a = zigzag(21)
    feat_b = zigzag(21,True)
    x,y = zip(*feat_b)
    ax.plot(x,y,color=c1,linewidth=0.5)
    x,y = zip(*feat_a)
    ax.plot(x,y,color=c2,linewidth=0.5)
    ax.set_title('worst case:\n orthogonal zigzags',y = -0.3)
    # display and save    
    
    
    folder = r"C:\CaGIS Board Dropbox\cantaloupe bob\Research\Projects\line difference metrics\Hausdorff\images"
    filename = "features_for_time_testing.png"
    plt.savefig(f"{folder}\\{filename}",dpi=300)
    plt.show()

# let's get the ball rolling with a simple test
vertex_counts = [20,50,100,200,500,1000,2000,5000,10000,20000]
vertex_counts = vertex_counts[:10]
brute_force = False



# create_example_plots()    
# test_same_polyline()    
time_test_unequal(cannonball(), vertex_counts, brute_force)
# time_test_zigzag(vertex_counts[:-4])