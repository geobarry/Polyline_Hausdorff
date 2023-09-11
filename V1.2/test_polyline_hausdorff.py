# -*- coding: utf-8 -*-
"""
Tests to verify accuracy of polyline Hausdorff distance calculation

Created on Mon May  3 10:14:29 2021

@author: bjkronenfeld
"""

import polyline_hausdorff as h
import utils_hausdorff as hu
from test_utils import case, plot_hausdorff_solution
import unittest
import random


class hausdorff_tests(unittest.TestCase):
    """Tests functions in the geom_utils.py module"""
    def setUp(self):
        self.verbose = False
        self.numits = 10;
    def test_transform_matches(self):
        # we are going to take all test cases, apply a randomly generated
        # affine transformation to both A & B,
        # and make sure that the Hausdorff distance between the transformed
        # polylines is the same as the original polylines
        def random_transform(A,B):
            # randomly shift, rescale and rotate points
            # shift
            dx = 10000 * (random.random()-0.5)
            dy = 10000 * (random.random()-0.5)
            A_prime = [(x+dx,y+dy) for x,y in A]
            B_prime = [(x+dx,y+dy) for x,y in B]
            return A_prime, B_prime
        # loop through all test cases
        print("working through {} random iterations".format(self.numits))
        for itnum in range(self.numits):
            print("iteration {}".format(itnum+1))
            for s in range(case(-1)):
                # obtain two polylines
                A,B = case(s)
                Aprime,Bprime = random_transform(A,B)
                # check the Hausdorff distance
                H,srcloc,srcline,srccomp,trglocs,trgcomps = h.polyline_hausdorff(A, B)
                Hprime,srclocprime,srclineprime,srccompprime,trglocsprime,trgcompsprime = h.polyline_hausdorff(Aprime, Bprime)
                self.assertAlmostEqual(H, Hprime)
                # check the Hausdorff distance from each segment
                def segHausdorff_wrapper(C,D,seg):
                    D_seg_idx = hu.seg_idx(D)
                    startNearSeg = hu.nearSegment(C, D, seg, D_seg_idx)
                    endNearSeg = hu.nearSegment(C, D, seg+1, D_seg_idx)
                    startcomp, startd = hu.nearComponent(C, D, seg, startNearSeg)
                    endcomp, endd = hu.nearComponent(C, D, seg+1, endNearSeg)
                    candidates = hu.candidateComponents(C,D,seg,startd,endd,D_seg_idx)
                    d,k,comps = h.hausdorff_segment(C, D, seg, startcomp, startd, endcomp, endd, candidates)
                    return d,k,comps
                def make_comparison(A,B,Aprime,Bprime,seg):
                    # original
                    d,k,comps = segHausdorff_wrapper(A,B,seg)
                    # transformed
                    dprime,kprime,compsprime = segHausdorff_wrapper(Aprime,Bprime,seg)
                    # compare
                    report = "Case: {} seg: {}".format(s,seg)
                    # visualize on error
                    if abs(d-dprime) > 0.001:
                        print("dprime: {}".format(dprime))
                        print("\nORIGINAL")
                        test_hausdorff_segment(A,B,seg,True,report + " original ")
                        print("\nTRANSFORMED")
                        test_hausdorff_segment(Aprime,Bprime,seg,True,report + " transformed ")
                    self.assertAlmostEqual(d,dprime,msg=report)
                # from A
                for seg in range(len(A)-1):
                    make_comparison(A,B,Aprime,Bprime,seg)
                # from B                
                for seg in range(len(B)-1):
                    make_comparison(B,A,Bprime,Aprime,seg)
            
def test_hausdorff_segment(A,B,a,verbose=False, titleprefix="",savefolder= ""):    
    B_idx = hu.seg_idx(B)
    startNearSeg = hu.nearSegment(A, B,a,B_idx)
    endNearSeg = hu.nearSegment(A, B, a+1,B_idx)
    startcomp, startd = hu.nearComponent(A, B, a, startNearSeg)
    endcomp, endd = hu.nearComponent(A, B, a+1, endNearSeg)
    candidates = hu.candidateComponents(A,B,a,startd,endd,B_idx)
    near_comps = h.segment_traversal(A, B, a, startcomp, startd, endcomp, endd, candidates,verbose)
    d,k,hausdorff_comps = h.hausdorff_segment(near_comps, verbose)
    msg = ""
    
    if True:
        # show sequence of nearest components
        msg += " -> ".join([hu.component_label(near_comps[i][2])for i in range(len(near_comps))])
    elif len(hausdorff_comps)==1:
        # show Hausdorff component
        msg += "to {}".format(hu.component_label(hausdorff_comps[0]))
    else:
        msg += "to {} and {}".format(hu.component_label(hausdorff_comps[0]),hu.component_label(hausdorff_comps[1]))
    msg += '   (dist: {:.3f})'.format(d)
    try:
        x = A[a][0] + k * (A[a+1][0]-A[a][0])
    except:
        print("k: {}".format(k))
    y = A[a][1] + k * (A[a+1][1]-A[a][1])
    imagefile = ""
    if savefolder != "" and savefolder != None:
        imagefile = r"{}\{}.png".format(savefolder,titleprefix)
    plot_hausdorff_solution(A, a, B, title = "{}:  {}".format(titleprefix, msg), imagefile = imagefile, show_labels = True)


def test_all_hausdorff_segments(A,B,verbose=False,titleprefix="",do_reverse = False,savefolder=None):
    for a in range(len(A)-1):
        print("polyline A segment {}...".format(a))
        test_hausdorff_segment(A,B,a,verbose,titleprefix + "-{}".format(a),savefolder)
    if do_reverse:
        for b in range(len(B)-1):
            print("polyline B segment {}...".format(b))
            test_hausdorff_segment(B,A,b,verbose,titleprefix + "-{}-rev".format(b),savefolder)

def printComponentInfo(A,B,a,comps):
    # print out switch points
    print("\n*** switch points ***")
    colheads = ""
    for c in comps:
        colheads += "\t{}".format(hu.component_label(c))
    t1msgs = []
    t2msgs = []
    for c1 in comps:
        dr1 = hu.distanceRepresentation(A, B, a, c1)
        t1msg = hu.component_label(c1)
        t2msg = hu.component_label(c1)
        for c2 in comps:
            dr2 = hu.distanceRepresentation(A, B, a, c2)
            swpts = hu.switchPoint(dr1, dr2)
            while len(swpts) < 2:
                swpts.append("-")
            t1msg += "\t{}".format(swpts[0])
            t2msg += "\t{}".format(swpts[1])
        t1msgs.append(t1msg)
        t2msgs.append(t2msg)
    print(colheads)
    for msg in t1msgs:
        print(msg)
    print("")
    print(colheads)
    for msg in t2msgs:
        print(msg)

    # print out effective intervals
    print("\n*** effective intervals ***")
    for c in comps:
        cei = hu.effectiveInterval(A, B, a, c)
        print("{}: {}".format(hu.component_label(c),cei))
    
    print("\n")

def test_all_cases_visually(savefolder = None):
    for s in range(case(-1)):
        print("CASE {}".format(s))
        A,B = case(s)
        test_all_hausdorff_segments(A, B, verbose=False, titleprefix="Case {}".format(s),savefolder = savefolder)
    print("\nTO CHECK RESULTS:")
    print("Review all diagrams. The red point should be the location on the")
    print("highlighted orange segment that is furthest from the blue polyline.")
    print("The gray circle should touch the blue polyline at one or two ")
    print("locations, indicating the location(s) on the blue polyline that are ")
    print("closest to the red point.")

def test_case_visually(case_number, seg, reverse = False):
    if reverse:
        B,A = case(case_number)
        direction = "_rev"
    else:
        A,B = case(case_number)
        direction = ""
    msg = "case {}{}_s{}".format(case_number,direction,seg)
    test_hausdorff_segment(A,B,seg,True,msg)

def test_overlap():
    A,B = case(2)
    verbose = True
    test_hausdorff_segment(A,B,1,verbose,"overlap ")



#unittest.main()
savefolder = r"C:\CaGIS Board Dropbox\cantaloupe bob\Research\Projects\line_difference_metrics\Hausdorff\images\test_cases_v2"
               
savefolder = r"C:\temp"               

# test_all_cases_visually(savefolder = savefolder)

A,B = case(6)

plot_hausdorff_solution(A, 0, B, show_labels = True, k_use_hausdorff=False,k=0.15)

# test_case_visually(0,0,False) # case from Hangouet
# test_case_visually(9,0)

# A,B = case(9)
# hu.distanceRepresentation(A, B, 0, (True,3))

# test_case_visually(0,0) # going backwards from vertex to segment
# test_case_visually(4,1,False) # navigating along the straight set of segments
# test_case_visually(6,2,True)

# test_case_visually(2,1)
# test_case_visually(2,2,True)

# test_case_visually(6,0,False) # situation where we have to double back
# test_case_visually(4,0,True)

# test_case_visually(5,1) # navigating along straight set of segments
# test_overlap()
print("\nfinished.")