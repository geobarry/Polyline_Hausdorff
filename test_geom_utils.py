# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
 Name:        geom utils
 Purpose:     various geometry functions
 License:     MIT License
 Updated:     April 2021
 Notes:       Docstrings follow numpy format 
              https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard
-------------------------------------------------------------------------------
"""

import unittest
import utils_geom as g
import math
import random

class geom_utils_tests(unittest.TestCase):
    """Tests functions in the geom_utils.py module"""
    def setUp(self):
        self.verbose = False

    def test_angle(self):
        # TEST ANGLE FUNCTION
        if self.verbose:
            print("************************")
            print("TESTING ANGLE FUNCTION")
        # test 1: series of angles radially around a fulcrum
        if self.verbose:
            print("---RADIAL SERIES---")
        fulcrum = (100,100)
        a = (100,200)
        b_list = [(100,150),(150,150),(200,100),(150,50),(100,0),(25,25),(-50,100),(0,200)]
        deg_list = [0,45,90,135,180,135,90,45] # known correct values
        for i in range(len(b_list)):
            b = b_list[i]
            theta = g.angle([fulcrum,a],[fulcrum, b])
            if self.verbose:
                print("Target: {}  Computed: {}".format(deg_list[i],math.degrees(theta)))
            self.assertAlmostEqual(deg_list[i],math.degrees(theta),msg="radial series")
            
        # test 2: series of random angles
        if self.verbose:
            print("---RANDOM ANGLES---")
        for i in range(10):
            # get random angles
            angle1 = random.random()*math.pi*2
            angle2 = random.random()*math.pi*2
            # answer is difference between them
            answer = abs(angle1-angle2)
            if answer > math.pi:
                answer = 2*math.pi-answer
            # construct vectors with random start
            v1a = (100*random.random(),100*random.random())
            dx = 100*math.cos(angle1)
            dy = 100*math.sin(angle1)
            v1b = (v1a[0]+dx,v1a[1]+dy)
            v1=[v1a,v1b]
            v2a = (100*random.random(),100*random.random())
            dx = 100*math.cos(angle2)
            dy = 100*math.sin(angle2)
            v2b = (v2a[0]+dx,v2a[1]+dy)
            v2=[v2a,v2b]
            func_out = g.angle(v1,v2)
            if self.verbose:
                print("Target: {}  Computed: {}".format(answer,func_out))
            self.assertAlmostEqual(answer,func_out,msg="random angles")

    def test_area(self):
        # TEST AREA FUNCTION
        if self.verbose:
            print("************************")
            print("TESTING AREA FUNCTION")
        # test 1: series of well-known geometries
        if self.verbose:
            print("--WELL KNOWN GEOMETRIES--")
        square = [(0,0),(0,100),(100,100),(100,0),(0,0)]
        A = g.area(square)
        if self.verbose:
            print("Square >> Target: {}  Computed {}".format(10000,A))
        self.assertAlmostEqual(10000,A,msg="Square")
        rectangle = [(-10,0),(-10,100),(990,100),(990,0),(-10,0)]
        A = g.area(rectangle)
        if self.verbose:
            print("Rectangle >> Target: {} Computed {}".format(100000,A))
        self.assertAlmostEqual(100000,A,msg="Rectangle")
        sqrt3 = math.sqrt(3)
        hexagon = [(0,1),(sqrt3,0.5),(sqrt3,-0.5),(0,-1),(-sqrt3,-0.5),(-sqrt3,0.5),(0,1)]
        A = g.area(hexagon)
        if self.verbose:
            print("Rectangle >> Target: {} Computed {}".format(sqrt3*3,A))
        self.assertAlmostEqual(sqrt3*3,A,msg="Rectangle")
        # test 2: series of random triangles
        if self.verbose:
            print("--RANDOM TRIANGLES--")
        for i in range(10):
            # get random points for two sides
            x1 = 100*random.random()
            y1 = 100*random.random()
            x2 = 100*random.random()
            y2 = 100*random.random()
            # get random distance (height) from base
            h = 100*random.random()
            # target area is 1/2 base x height
            dx = x2-x1
            dy = y2-y1
            base = math.sqrt(dx**2 + dy**2)
            target = base * h/2
            # get random location along line thru base
            k = -3 + random.random()*7
            xk = x1 + k*dx # x-coordinate of location k
            yk = y1 + k*dy # y-coordinate of location k
            # place third point at distance d perpendicular from location k
            # perpendicular direction is (dy,-dx)
            dx,dy = dy,-dx
            # scale to desired height
            dx = dx*h/base
            dy = dy*h/base
            # get third point
            x3 = xk + dx
            y3 = yk + dy
            # calculate area
            A = g.area([(x1,y1),(x2,y2),(x3,y3)])
            if self.verbose:
                print("Target: {} Computed {}".format(target,A))
            self.assertAlmostEqual(target,A,msg="random triangle{}".format(i+1))

    def test_distance(self):
        # TEST DISTANCE FUNCTION
        if self.verbose:
            print("************************")
            print("TESTING DISTANCE FUNCTION")
        # test 1: series of pythagorean triples
        if self.verbose:
            print("---PYTHAGOREAN TRIPLES---")
        pt1 = (0,0)
        pt2 = (3,4)
        if self.verbose:
            print("Target: {} Computed {}".format(5,g.distance(pt1, pt2)))
        self.assertAlmostEqual(5,g.distance(pt1,pt2))
        pt1 = (0,0)
        pt2 = (-5,-12)
        if self.verbose:
            print("Target: {} Computed {}".format(13,g.distance(pt1, pt2)))
        self.assertAlmostEqual(13,g.distance(pt1,pt2))
        xbase=random.random()
        ybase=random.random()
        pt1 = (xbase,ybase)
        pt2 = (xbase+8,ybase+15)
        if self.verbose:
            print("Target: {} Computed {}".format(17,g.distance(pt1, pt2)))
        self.assertAlmostEqual(17,g.distance(pt1,pt2))
        # test 2: random horizontal and vertical lines
        if self.verbose:
            print("---RANDOM HORIZONTAL & VERTICAL LINES---")
        # horizontal
        for i in range(5):
            xbase=random.random()*100
            ybase=random.random()*100
            d = random.random()*100
            pt1 = (xbase,ybase)
            pt2 = (xbase, ybase + d)
            pt3 = (xbase + d, ybase)
            if self.verbose:
                print("Target: {} Computed {}".format(d,g.distance(pt1, pt2)))
                print("Target: {} Computed {}".format(d,g.distance(pt1, pt3)))
            self.assertAlmostEqual(d,g.distance(pt1,pt2))
            self.assertAlmostEqual(d,g.distance(pt1,pt3))

    def test_distance_to_line(self):
        if self.verbose:
            print("************************")
            print("TESTING DISTANCE_TO_LINE FUNCTION")
        # test 1: horizontal and vertical cases
        if self.verbose:
            print("--HORIZONTAL LINE--")
        a = (0,0)
        b = (100,0)
        p = (50,500)
        if self.verbose:
            print("Target: {} Computed {}".format(500,g.distance_to_line(p,a,b)))
        self.assertAlmostEqual(500,g.distance_to_line(p,a,b))
        p = (150,-500)
        if self.verbose:
            print("Target: {} Computed {}".format(500,g.distance_to_line(p,a,b)))
        self.assertAlmostEqual(500,g.distance_to_line(p,a,b))
        if self.verbose:
            print("--VERTICAL LINE--")
        a = (0,0)
        b = (0,100)
        p = (500,50)
        if self.verbose:
            print("Target: {} Computed {}".format(500,g.distance_to_line(p,a,b)))
        self.assertAlmostEqual(500,g.distance_to_line(p,a,b))
        p = (-500,50)
        if self.verbose:
            print("Target: {} Computed {}".format(500,g.distance_to_line(p,a,b)))
        self.assertAlmostEqual(500,g.distance_to_line(p,a,b))
        # test 2: random offsets
        if self.verbose:
            print("--RANDOM OFFSETS--")
        for i in range(10):
            # get random points for two sides
            x1 = 100*random.random()
            y1 = 100*random.random()
            x2 = 100*random.random()
            y2 = 100*random.random()
            # get random distance (height) from base
            h = 100*random.random()
            # target area is 1/2 base x height
            dx = x2-x1
            dy = y2-y1
            base = math.sqrt(dx**2 + dy**2)
            # get random location along line thru base
            k = -3 + random.random()*7
            xk = x1 + k*dx # x-coordinate of location k
            yk = y1 + k*dy # y-coordinate of location k
            # place third point at distance d perpendicular from location k
            # perpendicular direction is (dy,-dx)
            dx,dy = dy,-dx
            # scale to desired height
            dx = dx*h/base
            dy = dy*h/base
            # get third point
            x3 = xk + dx
            y3 = yk + dy
            # test
            a = (x1,y1)
            b = (x2,y2)
            p = (x3,y3)
            d = g.distance_to_line(p,a,b)
            if self.verbose:
                print("Target: {} Computed {}".format(h,d))
            self.assertAlmostEqual(h,d)

    def test_intersection(self):
        if self.verbose:
            print("************************")
            print("TESTING INTERSECTION FUNCTION")
        # test 1: horizontal and vertical lines
        if self.verbose:
            print("--HORIZONTAL AND VERTICAL LINE INTERSECTIONS--")
        a = (0,0)
        b = (100,0)
        c = (50,-50)
        d = (50,50)
        p = []
        p.append(g.intersection(a,b,c,d))
        p.append(g.intersection(b,a,c,d))
        p.append(g.intersection(a,b,d,c))
        p.append(g.intersection(b,a,d,c))
        for q in p:
            if self.verbose:
                print("Target: {} Computed {}".format((50,0),q))
            self.assertAlmostEqual(50,q[0])
            self.assertAlmostEqual(0,q[1])
        
        # test 2: parallel lines
        if self.verbose:
            print("--PARALLEL LINES--")
        a = (0,0)
        b = (100,0)
        c = (50,-50)
        d = (150,-50)
        if self.verbose:
            print("Target: {} Computed {}".format("*No Intersection*",g.intersection(b,a,d,c)))
        self.assertEqual(g.intersection(b,a,d,c),(None,None))
        a = (0,0)
        b = (100,50)
        c = (50,-50)
        d = (150,0)
        if self.verbose:
            print("Target: {} Computed {}".format("*No Intersection*",g.intersection(b,a,d,c)))
        self.assertEqual(g.intersection(b,a,d,c),(None,None))
    
        # test 3: random lines through a point
        if self.verbose:
            print("--RANDOM LINES THROUGH A POINT--")
        for i in range(10):
            # get three random points
            crosspt = (1000*random.random(),1000*random.random())
            a =  (1000*random.random(),1000*random.random())
            b =  (1000*random.random(),1000*random.random())
            # use k-values to get second points on each line
            # get two random "k-values"
            ka = -3 + 7*random.random()
            kb = -3 + 7*random.random()
            # get dx and dy for each line
            dxa = a[0]-crosspt[0]
            dya = a[1]-crosspt[1]
            dxb = b[0]-crosspt[0]
            dyb = b[1]-crosspt[1]
            # move along line from crosspt to new point using k-value
            a2 = (crosspt[0]+ka*dxa,crosspt[1]+ka*dya)
            b2 = (crosspt[0]+kb*dxb,crosspt[1]+kb*dyb)
            p = g.intersection(a,a2,b,b2)
            if self.verbose:
                print("Target: {} Computed {}".format(crosspt,p))
            self.assertAlmostEqual(crosspt[0],p[0])
            self.assertAlmostEqual(crosspt[1],p[1])
    
    def test_is_monotonic(self):
        test_vals = []
        test_vals.append(((1,2,3),True))
        test_vals.append(((1,3,2),False))
        test_vals.append(((1,2,2),True))
        test_vals.append(((2,1,3),False))
        test_vals.append(((2,3,1),False))
        test_vals.append(((2,1,1),True))
        test_vals.append(((2,3,3),True))
        test_vals.append(((2,2,1),True))
        test_vals.append(((2,2,3),True))
        test_vals.append(((3,1,2),False))
        test_vals.append(((3,2,1),True))
        test_vals.append(((3,2,2),True))
        for tval in test_vals:
            a,b,c = tval[0]
            result = g.is_monotonic(a, b, c)
            target = tval[1]
            if self.verbose:
                print("Target: {} ... Result: {}".format(target,result))
            self.assertAlmostEqual(target,result)
    
    def test_project_pt_to_line(self):
        if self.verbose:
            print("************************")
            print("TESTING PROJECT_PT_TO_LINE FUNCTION")
        # test 1: horizontal and vertical lines
        if self.verbose:
            print("--HORIZONTAL LINE--")
        a = (0,0)
        b = (100,0)
        p = (50,500)
        prj = g.project_pt_to_line(p,a,b)
        
        if self.verbose:
            print("Target: {} Computed {}".format((50,0),prj))
        self.assertAlmostEqual(50,prj[0])
        self.assertAlmostEqual(0,prj[1])
    
        p = (150,-500)
        prj = g.project_pt_to_line(p,a,b)

        if self.verbose:
            print("Target: {} Computed {}".format((150,0),prj))
        self.assertAlmostEqual(150,prj[0])
        self.assertAlmostEqual(0,prj[1])

        if self.verbose:
            print("--VERTICAL LINE--")
        a = (0,0)
        b = (0,100)
        p = (500,50)

        prj = g.project_pt_to_line(p,a,b)
        if self.verbose:
            print("Target: {} Computed {}".format((0,50),prj))
        self.assertAlmostEqual(0,prj[0])
        self.assertAlmostEqual(50,prj[1])

        p = (-500,50)
        prj = g.project_pt_to_line(p,a,b)
        if self.verbose:
            print("Target: {} Computed {}".format((0,50),prj))
        self.assertAlmostEqual(0,prj[0])
        self.assertAlmostEqual(50,prj[1])

        # test 2: random offsets
        if self.verbose:
            print("--RANDOM OFFSETS--")
        for i in range(10):
            # get random points for two sides
            x1 = 100*random.random()
            y1 = 100*random.random()
            x2 = 100*random.random()
            y2 = 100*random.random()
            # get random distance (height) from base
            h = 100*random.random()
            # target area is 1/2 base x height
            dx = x2-x1
            dy = y2-y1
            base = math.sqrt(dx**2 + dy**2)
            # get random location along line thru base
            k = -3 + random.random()*7
            xk = x1 + k*dx # x-coordinate of location k
            yk = y1 + k*dy # y-coordinate of location k
            # place third point at distance d perpendicular from location k
            # perpendicular direction is (dy,-dx)
            dx,dy = dy,-dx
            # scale to desired height
            dx = dx*h/base
            dy = dy*h/base
            # get third (unprojected) point
            x3 = xk + dx
            y3 = yk + dy
            # test
            a = (x1,y1)
            b = (x2,y2)
            p = (x3,y3)
            prj = g.project_pt_to_line(p,a,b)
            trg = (xk,yk)
            if self.verbose:
                print("Target: {} Computed {}".format(trg,prj))
            self.assertAlmostEqual(trg[0],prj[0])
            self.assertAlmostEqual(trg[1],prj[1])


    
    def test_kvalue(self):
        if self.verbose:
            print("************************")
            print("TESTING PROJECT_PT_TO_LINE FUNCTION")
            print("--NORMAL CASE--")
        a=(0,0)
        b=(10,10)
        p=(5,5)
        if self.verbose:
            print("Target: {} --- Computed: {}".format(0.5,g.kvalue(p, a, b)))
        self.assertAlmostEqual(0.5,g.kvalue(p,a,b))

        p=(-5,-5)
        if self.verbose:
            print("Target: {} --- Computed: {}".format(-0.5,g.kvalue(p, a, b)))
        self.assertAlmostEqual(-0.5,g.kvalue(p,a,b))


        p=(20,20)
        if self.verbose:
            print("Target: {} --- Computed: {}".format(2,g.kvalue(p, a, b)))
        self.assertAlmostEqual(2,g.kvalue(p,a,b))

        if self.verbose:
            print("--VERTICAL LINE SEGMENT--")
        a=(0,0)
        b=(0,10)
        p=(0,5)
        if self.verbose:
            print("Target: {} --- Computed: {}".format(0.5,g.kvalue(p, a, b)))
        self.assertAlmostEqual(0.5,g.kvalue(p,a,b))

        p=(0,-5)
        if self.verbose:
            print("Target: {} --- Computed: {}".format(-0.5,g.kvalue(p, a, b)))
        self.assertAlmostEqual(-0.5,g.kvalue(p,a,b))

        p=(0,20)
        if self.verbose:
            print("Target: {} --- Computed: {}".format(2,g.kvalue(p, a, b)))
            print("--HORIZONTAL LINE SEGMENT--")
        self.assertAlmostEqual(2,g.kvalue(p,a,b))

        a=(0,0)
        b=(10,0)
        p=(5,0)
        if self.verbose:
            print("Target: {} --- Computed: {}".format(0.5,g.kvalue(p, a, b)))
        self.assertAlmostEqual(0.5,g.kvalue(p,a,b))

        p=(-5,0)
        if self.verbose:
            print("Target: {} --- Computed: {}".format(-0.5,g.kvalue(p, a, b)))
        self.assertAlmostEqual(-0.5,g.kvalue(p,a,b))

        p=(20,0)
        if self.verbose:
            print("Target: {} --- Computed: {}".format(2,g.kvalue(p, a, b)))
            print("--SLIGHT FLOATING POINT ERRORS--")
        self.assertAlmostEqual(2,g.kvalue(p,a,b))

        a=(0.00000000000003,-0.00000000001)
        b=(10,10)
        p=(5.000000000001,4.99999999)
        if self.verbose:
            print("Target: {} --- Computed: {}".format(0.5,g.kvalue(p, a, b)))
        self.assertAlmostEqual(0.5,g.kvalue(p,a,b))
        
        p=(-4.9999999998,-5.000000002)
        if self.verbose:
            print("Target: {} --- Computed: {}".format(-0.5,g.kvalue(p, a, b)))
        self.assertAlmostEqual(-0.5,g.kvalue(p,a,b))

        p=(20.0000000002,19.99999999991)
        if self.verbose:
            print("Target: {} --- Computed: {}".format(2,g.kvalue(p, a, b)))
        self.assertAlmostEqual(2,g.kvalue(p,a,b))

unittest.main()