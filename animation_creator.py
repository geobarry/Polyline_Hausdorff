# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 14:28:36 2023

@author: Administrator
"""
import test_utils as tu

# run a test case and show in interpreter
A,B = tu.case(0)
tu.plot_hausdorff_solution(A, 1, B, 
                            show_labels = True, 
                            k=0.7)

# run a test case and save to a file
image_folder = r'C:\CaGIS Board Dropbox\cantaloupe bob\Research\Conferences\GIScience\GIScience2023\Cartography Workshop\images'
A,B = tu.case(6)
image_file = '{}\\{}.png'.format(image_folder,'case 6')
tu.plot_hausdorff_solution(A, 0, B, 
                            show_labels = True, 
                            k=0.75,
                            k_buffer = 0.08,
                            imagefile = image_file,
                            dpi = 400)
