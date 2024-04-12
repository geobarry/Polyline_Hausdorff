# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 15:37:16 2023

@author: Barry Kronenfeld 
"""

import integrals
import utils_hausdorff 
import math
def segment_integral_test_cases(id):
    """Returns parameters (L,k0,k1,dr) for integral function
        Parameters
        ----------
        id : integer
            if -1 then will return the total number of cases

    """
    cases = []
    # cosine thirty degrees
    cos30 = math.cos(math.radians(30))
    # parallel lines
    cases.append((5,0.2,0.6,(True,None,0.6))) # Answer should be 6
    # B segment does not span intersection
    cases.append((5,0.4,0.8,(True,1.2,cos30))) # Answer should be 3
    # B segment spans intersection
    cases.append((5,0.2,1.0,(True,0.6,cos30))) # Answer should be 2
    if id == -1:
        return len(cases)
    else:
        return cases[id]

def vertex_integral_test_cases(id):
    """Returns parameters (L,k0,k1,dr) for integral function
        Parameters
        ----------
        id : integer
            if -1 then will return the total number of cases

    """
    # (L,k0,k1,kv,q)
    cases = []
    # unit line with vertex at distance one in the middle
    cases.append((1,0,1,(False,0.5,1)))
    # unit line with vertex at distance zero in the middle
    cases.append((1,0,1,(False,0.5,0)))
    # unit line with vertex at distance zero off center
    cases.append((1,0.05,0.75,(False,0.5,0)))
    if id == -1:
        return len(cases)
    else:
        return cases[id]

def integral_brute_force(L,k0,k1,dr,n = 100):
    """
        Calculates the integral using the bruteforce method. 
        Arguments are the same as integrals.vertex_integral
    """
    totD = 0
    for i in range(n):
        k = k0 + ((i + 0.5)/n) * (k1 - k0)
        d = utils_hausdorff.componentDistance(dr, k, L)
        totD = totD + d
    avgD = totD/n
    area = avgD * abs(k1-k0) * L
    return area

def test_all_cases():
    for id in range(vertex_integral_test_cases(-1)):
        L,k0,k1,dr = vertex_integral_test_cases(id)        
        A = integrals.component_integral(L, k0, k1, dr)
        bruteforce_area = integral_brute_force(L, k0, k1, dr, 700)
        error = abs(bruteforce_area - A)
        print("case {}: {:.5f} (error: {:.6f})".format(id,A,error))
    for id in range(segment_integral_test_cases(-1)):
        L,k0,k1,dr = segment_integral_test_cases(id)
        A = integrals.component_integral(L, k0, k1, dr)
        bruteforce_area = integral_brute_force(L,k0,k1,dr, 10)
        error = abs(bruteforce_area - A)
        print("case {}: {:.5f} (error: {:.6f})".format(id,A,error))
            

test_all_cases()
