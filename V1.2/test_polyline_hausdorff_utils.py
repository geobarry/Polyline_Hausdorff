# -*- coding: utf-8 -*-
"""
Created on Wed May 26 13:30:22 2021

@author: bjkronenfeld
"""
import unittest
import matplotlib.pyplot as plt
import math as m
import utils_hausdorff as hu
import utils_geom as g

def drawLines(lineList,title="", lastFigure = False):
    """Plots a set of lines. To plot a point, input a line with one point."""
    for line in lineList:
        if len(line) > 1:
            x,y = zip(*line)
            plt.plot(x,y)
        else:
            plt.plot(line[0][0],line[0][1],'ro')
    plt.suptitle(title)
    if lastFigure == False:
        plt.figure()
    plt.show()

class Hausdorff_utils_tests(unittest.TestCase):
    """Tests the functions in Hausdorff_utils.py"""
    def setUp(self):
        self.verbose = False
        self.showPlots = False

    def test_componentDistance(self):
        print(">>> FUNCTION componentDistance")

        def report(answer,target,case):
            if self.verbose:
                print("\n{}".format(case))
                print("computed: {}".format(answer))
                print("target:   {}".format(target))
            self.assertAlmostEqual(answer,target,msg=case)

        A = [(100,150),(150,150),(200,100),(150,50),(100,0),(25,25),(-50,100),(0,200)]
        case = 'test case for vert Distance'
        dr = (False, 0.03, 1)
        k = 0.5
        a = 0
        # answer should be sqrt(0.47**2 + 1**2) * 50
        target = 50*m.sqrt(0.47**2 + 1)
        len_a = g.distance(A[a],A[a+1])
        answer = hu.componentDistance(dr, k, len_a)
        report(answer,target,case)
        
        case = 'test case for seg Distance'
        dr = (True, 0.03, 1)
        k = 0.5 
        a = 2
        # segment is vertica. answer should be 0.47 * length of segment
        target = 0.47 * m.sqrt(50**2+50**2)
        len_a = g.distance(A[a],A[a+1])
        answer = hu.componentDistance(dr, k, len_a)
        report(answer,target,case)
        if self.verbose:
            print("\n")

    def test_distanceRepresentation(self):
        def report(answer,target,case):
            if self.verbose:
                print("\n{}".format(case))
                print("computed: {}".format(answer))
                print("target:   {}".format(target))
            self.assertEqual(answer[0],target[0],msg=case + " - item 1")
            if answer[1] == None or target[1] == None:
                self.assertEqual(answer[1],target[1],msg=case + " - item 2")
            else:
                self.assertAlmostEqual(answer[1],target[1],msg=case + " - item 2")
            self.assertAlmostEqual(answer[2],target[2],msg=case + " - item 3")

        print(">>> FUNCTION distanceRepresentation")
        case="I am running the Vertex function"
        a = 1
        bcomp = (False, 2)
        A = [(100,14),(0,0),(1,0),(42,-1323123)]
        B = [(43,-.2),(-0.3,0.4),(86,42)]
        answer = hu.distanceRepresentation(A, B, a, bcomp)
        # result should be vertex 2 k & q
        target = (False,86,42)
        report(answer,target,case)
        
        case = "I am running the segment function"
        a = 3
        bcomp = (True,3)
        A = [(6,0),(3,4),(2,1),(0,0),(1,0)]
        B = [(7,1),(8,6),(5,6),(.6,.4),(1.3,.4)]
        answer = hu.distanceRepresentation(A, B, a, bcomp)
        # result should be segment 3 intersection k-value and sin theta
        target = (True,None,0.4)
        report(answer,target,case)
        if self.verbose:
            print("\n")

    def test_effectiveInterval(self):
        print(">>> FUNCTION effectiveInterval")
        def report(answer,target,case):
            if self.verbose:
                print("\n{}".format(case))
                print("computed: {}".format(answer))
                print("target:   {}".format(target))
            self.assertEqual(len(answer), len(target), msg=case + " - lenths unequal")
            if len(answer)==len(target):
                for i in range(len(answer)):
                    self.assertAlmostEqual(answer[i], target[i], msg = case + " - item i")

        case = "Input is a segment:"
        A = [(0,0), (10,0)]
        B = [(2,1), (3,1),(4,1)]
        a = 0
        bcomp = (True, 1)
        test_answer = hu.effectiveInterval(A, B, a, bcomp)
        target = (0.3,0.4)
        report(test_answer,target,case)
    
        case = "Input is a vertex:"
        bcomp = (False, 1)
        test_answer = hu.effectiveInterval(A, B, a, bcomp)
        target = (0.3,0.3)
        report(test_answer,target,case)

        
    def test_segDistance(self):
        print(">>> FUNCTION segDistance")
        def report_case(dr,k,len_a,target,case):
            computed = hu.segDistance(dr,k,len_a)
            if self.verbose:
                print('\n' + case)
                print("computed: {}\ntarget: {}".format(computed, target))
            self.assertAlmostEqual(computed, target)
                
        case = "--- vertical segment ---"
        dr = (True, 0.2, 1)
        k = 0.5 
        len_a = 2.5
        correct_answer = 0.75
        report_case(dr,k,len_a,correct_answer,case)
    
        case = "--- 45 degree angle ---"
        dr = (True, 0, 1/(2**0.5))
        k = 0.5 
        len_a = 100
        correct_answer = 50/(2**0.5)
        report_case(dr,k,len_a,correct_answer,case)
        
        case = "--- parallel segments ---"
        dr = (True, None, 0.5)
        k = 0.5 
        len_a = 100
        correct_answer = 50
        report_case(dr,k,len_a,correct_answer,case)
        if self.verbose:
            print("\n")        
    
    def test_segDistRep(self):
        print(">>> FUNCTION segDistRep")
        # Test for no intersection
        a = 0
        b = 0
        A = [(0,0),(2,0)]
        B = [(.6,.4),(1.3,.4)]
        answer = hu.segDistRep(A, B, a, b)
        target = (True,None,0.2)
        if self.verbose:
            print("\nCase 1: Lines do not intersect.")
            print("computed: {}".format(answer))
            print("target:   {}".format(target))
        if self.showPlots:
            drawLines([A, B], "segDistRep case 1: Lines do not intersect\n{}".format(answer))
        [self.assertAlmostEqual(answer[i],target[i], msg="case 1 value {}".format(i)) for i in range(3)]
        # Test for an intersection
        B[0] = (0.6,0.8)
        answer = hu.segDistRep(A, B, a, b)
        angle = m.atan(0.4/0.7)
        target = (True,1,m.sin(angle))
        if self.verbose:
            print("\nCase 2: Lines intersect.")
            print("computed: {}".format(answer))
            print("target:   {}".format(target))
        if self.showPlots:
            drawLines([A, B], "segDistRep case 2: Lines intersect\n{}".format(answer))
        [self.assertAlmostEqual(answer[i],target[i], msg="case 2 value {}".format(i)) for i in range(3)]
        if self.verbose:
            print("\n")

    def test_segEffectiveInterval(self):
        print(">>> FUNCTION segEffectiveInterval")

        a = 0
        b = 0
        A = [(0,0),(1,0)]
        B = [(0.4,0.5),(0.9,0.2)]
        t= hu.segEffectiveInterval(A,B,a,b)
        if self.verbose:
            print("\noriginal test case:")
            print("computed: {}".format(t))
            #print("target:   {}".format(target))
    
        msg = "b slopes downward 45 degrees"
        A = [(0,0),(10,0)]
        B = [(4,3),(6,1)]
        a = 0 # first segment of A
        b = 0 # first segment of B
        correct_answer = (0.1,0.5)
        ei = hu.segEffectiveInterval(A,B,a,b)
        if self.verbose:
            print("\n" + msg)
            print("computed: {}".format(ei))
            print("target:   {}".format(correct_answer))
        [self.assertAlmostEqual(correct_answer[i],ei[i], msg=msg) for i in range(2)]
    
        A = [(0,0),(10,0)]
        B = [(4,1),(6,3)]
        a = 0 # first segment of A
        b = 0 # first segment of B
        correct_answer = (0.5,0.9)
        ei = hu.segEffectiveInterval(A,B,a,b)
        if self.verbose:
            print("\nb slopes upward 45 degrees")
            print("computed: {}".format(ei))
            print("target:   {}".format(correct_answer))
        [self.assertAlmostEqual(correct_answer[i],ei[i], msg="b slopes upward 45 degrees") for i in range(2)]
    
        A = [(0,0),(10,0)]
        B = [(4,1),(6,1)]
        a = 0 # first segment of A
        b = 0 # first segment of B
        correct_answer = (0.4,0.6)
        ei = hu.segEffectiveInterval(A,B,a,b)
        if self.verbose:
            print("\nb parallel to a")
            print("computed: {}".format(ei))
            print("target:   {}".format(correct_answer))
        [self.assertAlmostEqual(correct_answer[i],ei[i], msg="b parallel to a") for i in range(2)]
    
        A = [(0,0),(10,10)]
        B = [(3,9),(5,9)]
        a = 0 # first segment of A
        b = 0 # first segment of B
        correct_answer = (0.3,0.5)
        ei = hu.segEffectiveInterval(A,B,a,b)
        if self.verbose:
            print("\na rotated, b slopes downward 45 degrees")
            print("computed: {}".format(ei))
            print("target:   {}".format(correct_answer))
        [self.assertAlmostEqual(correct_answer[i],ei[i], msg="ba rotated, b slopes downward 45 degrees") for i in range(2)]
        if self.verbose:
            print("\n")

    def switchPoint_test(self):
        print(">>> FUNCTION vertDistance")
        def report_case(dr1,dr2,target,case):
            computed = hu.switchPoint(dr1, dr2)
            if self.verbose:
                print('\n' + case)
                print("computed: {}\ntarget: {}".format(computed, target))
            self.assertEqual(len(computed),len(target))
            if len(computed) == len(target):
                for i in range(len(target)):    
                    self.assertAlmostEqual(computed[i], target[i],2)
        
        case = "Vertex-Segment"
        dr1 = (False, 0.2, 0.1)
        dr2 = (True, 0.8, 0.5)
        target = [-0.383,0.383]
        report_case(dr1,dr2,target,case)
        
        case = "Segment-Segment"
        dr1 = (True, 2.5, 1)
        dr2 = (True, 2, 1)
        target = [2.25]
        report_case(dr1,dr2,target,case)
        
        case = "Vertex-Vertex"
        dr1 = (False, -0.2, 0.3)
        dr2 = (False, 0.7, 0.6)
        target=[0.4]
        report_case(dr1,dr2,target,case)
        
        case = "Segment-Vertex"
        dr2 = (False, 0.2, 0.1)
        dr1 = (True, 0.8, 0.5)
        target = [-0.383,0.383]
        report_case(dr1,dr2,target,case)
        
        if self.verbose:
            print("\n")
            
    def test_vertDist(self):
        print(">>> FUNCTION vertDistance")
        dr_list = [(True, 0.03, 1), (True, 0.05,0.5), (True, 0.09, 1.7), (True, 1.2, 2)]
        k_list = [0.5, 1, 1.2, 0.7]
        len_a_list = [2.5, 3, 0.07, 1]
        correct_answers = [2.76236,3.22064,0.14212,2.06155] # from MS Excel
        for i in range(len(k_list)):
            test_answer = hu.vertDistance(dr_list[i], k_list[i], len_a_list[i]) 
            if self.verbose:
                print("Target: {}\n Computed: {}".format(correct_answers[i],test_answer))
            self.assertAlmostEqual(correct_answers[i], test_answer,4)
        if self.verbose:
            print("\n")
            
    def test_vertdistrep(self):
        print(">>> FUNCTION vertDistRep")
        a = 0
        b = 0
        A = [(0,0),(3,0),(2,-1)]
        B = [(-0.3,0.4),(86,42)]
        answer = hu.vertDistRep(A,B,a,b)
        dA = g.distance(A[a], A[a+1])
        target = (False,-0.3/dA,0.4/dA)
        if self.verbose:
            print("computed: {}".format(answer))
            print("target:   {}".format(target))
        if self.showPlots:
            drawLines([A,B[:1]],"vertDistRep",lastFigure = True)
        [self.assertAlmostEqual(answer[i],target[i], msg="case 2 value {}".format(i)) for i in range(3)]
        if self.verbose:
            print("\n")

    def test_vertEffectiveInterval(self):
        print(">>> FUNCTION vertEffectiveInterval")

        def test_case(answer,target,case):
            if self.verbose:
                print("\n{}".format(case))
                print("computed: {}".format(answer))
                print("target:   {}".format(target))
            self.assertEqual(len(target),len(answer))
            [self.assertAlmostEqual(answer[i],target[i], msg=case) for i in range(len(target))]        
        
        case = 'Case 1a: start vertex left side'        
        a = 0
        b = 0
        A = [(0,0),(10,0)]
        B = [(6,4),(8,2)]
        target = (0,0.2)
        answer = hu.vertEffectiveInterval(A, B, a, b)
        test_case(answer,target,case)

        case = 'Case 1b: start vertex right side'        
        a = 0
        b = 0
        A = [(0,0),(10,0)]
        B = [(6,2),(4,2)]
        target = (0.6,1)
        answer = hu.vertEffectiveInterval(A, B, a, b)
        test_case(answer,target,case)

        case = 'Case 2: end vertex right side'        
        a = 0
        b = 1
        A = [(0,0),(10,0)]
        B = [(6,2),(4,2)]
        target = (0,0.4)
        answer = hu.vertEffectiveInterval(A, B, a, b)
        test_case(answer,target,case)

        case = 'Case 3: nbrs on left'        
        a = 0
        b = 1
        A = [(0,0),(10,0)]
        B = [(4,2),(6,3),(4,5)]
        target = (0.75,1)
        answer = hu.vertEffectiveInterval(A, B, a, b)
        test_case(answer,target,case)

        case = 'Case 4: nbrs on right'        
        a = 0
        b = 1
        A = [(0,0),(10,0)]
        B = [(6,2),(4,3),(6,5)]
        target = (0,0.25)
        answer = hu.vertEffectiveInterval(A, B, a, b)
        test_case(answer,target,case)

        case = 'Case 5: downward pointing vertex'        
        a = 0
        b = 1
        A = [(0,0),(10,0)]
        B = [(3,4),(5,2),(7,4)]
        target = (0.3,0.7)
        answer = hu.vertEffectiveInterval(A, B, a, b)
        test_case(answer,target,case)

        case = 'Case 6: upward pointing vertex'        
        a = 0
        b = 1
        A = [(0,0),(10,0)]
        B = [(3,4),(5,7),(7,4)]
        target = (float('-inf'),float('-inf'))
        answer = hu.vertEffectiveInterval(A, B, a, b)
        test_case(answer,target,case)
        if self.verbose:
            print("\n")
            
    def test_segSegSwitchPoint(self):
        print(">>> FUNCTION segSegSwitchPoint")
        def report_case(dr1,dr2,target,case):
            computed = hu.segSegSwitchPoint(dr1, dr2)
            if self.verbose:
                print('\n' + case)
                print("computed: {}\ntarget: {}".format(computed, target))
            self.assertEqual(len(computed),len(target))
            for i in range(len(target)):    
                self.assertAlmostEqual(computed[i], target[i],2)

        case = "Test1: Second segment of b parallel to a"
        dr1 = (True, -0.4, 0.476)
        dr2 = (True, None, 0.4141)
        target = [-1.26996,0.469958]
        report_case(dr1,dr2,target,case)
        
        case = "Test2: Both segments of b parallel to a"
        dr1 = (True, None, 0.476)
        dr2 = (True, None, 0.4141)
        target = []
        report_case(dr1,dr2,target,case)
        
        case = "Test3: First segment of b parallel to a"
        dr1 = (True, None, 0.476)
        dr2 = (True, -0.4, 0.4141)
        target = [-1.54948,0.749481]
        report_case(dr1,dr2,target,case)
        
        case = "Test4/Case3: Both segments of b make same angle with a"
        dr1 = (True, 2.5, 1)
        dr2 = (True, 2, 1)
        target = [2.25]
        report_case(dr1,dr2,target,case)
        
        case = "Test 5/Case 4: Normal case"
        dr1 = (True, -0.4, 0.476)
        dr2 = (True, 0.5, 0.9735)
        target = [0.20445,1.361106]
        report_case(dr1,dr2,target,case)
        if self.verbose:
            print("\n")
            
    def test_vertSegSwitchPoint(self):
        print(">>> FUNCTION vertSegSwitchPoint")
        def report_case(vdr,sdr,target,case):
            computed = hu.vertSegSwitchPoint(vdr, sdr)
            if self.verbose:
                print('\n' + case)
                print("computed: {}\ntarget: {}".format(computed, target))
            self.assertEqual(len(computed),len(target))
            for i in range(len(target)):    
                self.assertAlmostEqual(computed[i], target[i],2)
        
        # normal case
        case = "Normal case (from guide)"
        vdr = (False, 0.2, 0.1)
        sdr = (True, 0.8, 0.5)
        target = [-0.383,0.383]
        report_case(vdr,sdr,target,case)
        
        # no switch point
        case = "Case where segment is always closer than vertex"
        vdr = (False, 0.3, 0.5)
        sdr = (True, 0.5, 1/(2**0.5))
        target = [] 
        report_case(vdr,sdr,target,case)
        
        # vertex on segment; should return one point
        case = "Vertex is on segment"
        vdr = (False, 0.55, 0.6)
        sdr = (True, 1.35, 0.6)
        target = [0.1] 
        report_case(vdr,sdr,target,case)
    
        # perpendicular segment
        case = "Case of perpendicular segment"
        vdr = (False, 0.2, 0.3)
        sdr = (True, 1.1, 1)
        target = [0.6]
        report_case(vdr,sdr,target,case)
        
        # parallel segment
        case = "Case of parallel segment"
        vdr = (False, 0.1, 0.3)
        sdr = (True, None, 0.5)
        target = [-0.3,0.5]
        report_case(vdr,sdr,target,case)
        
        if self.verbose:
            print("\n")

    def test_vertVertSwitchPoint (self):
        print(">>> FUNCTION vertVertSwitchPoint")
        def report_case(dr1,dr2,target,case):
            computed = hu.vertVertSwitchPoint(dr1, dr2)
            if self.verbose:
                print('\n' + case)
                print("computed: {}\ntarget: {}".format(computed, target))
            self.assertAlmostEqual(computed, target)

        case = 'test values k and q are different'
        dr1 = (False, -0.2, 0.3)
        dr2 = (False, 0.7, 0.6)
        target=[0.4]
        report_case(dr1,dr2,target,case)
        
        case = 'test values for k are different, but q is equal'
        dr1 = (False, -0.2, 0.6)
        dr2 = (False, 0.7, 0.6)
        target=[0.25]
        report_case(dr1,dr2,target,case)

        case = 'test values for k are equal, but q is different'
        dr1 = (False, 0.7, 0.3)
        dr2 = (False, 0.7, 0.6)
        target=[]
        report_case(dr1,dr2,target,case)
        
        if self.verbose:
            print("\n")
        
unittest.main()