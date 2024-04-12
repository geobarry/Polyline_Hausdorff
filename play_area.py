# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 16:28:23 2023

@author: Administrator
"""

import test_utils as tu
A,B = tu.case(8)
seg_num = 0
tu.plot_hausdorff_solution(A, seg_num, B, 
                    show_labels = True, 
                    k=0.87,
                    k_buffer = 0,
                    dpi = 150)
