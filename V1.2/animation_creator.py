# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 14:28:36 2023

@author: Adam Camerer
"""
import test_utils as tu
import multiprocessing as mp
import sys
import os
import numpy as np
import time
import imageio
import shutil

# run a test case and show in interpreter
#A,B = tu.case(0)
#tu.plot_hausdorff_solution(A, 1, B, 
#                            show_labels = True, 
#                            k=0.7)

# run a test case and save to a file

#set folder and base file path
case_nums = [0, 4, 6, 7, 9]
base_image_folder = '/home/acamerer/Desktop/My_Project/Workflows/hausdroff_animation/outputs'

#wipe outputs folder
if os.path.isdir(base_image_folder):
        shutil.rmtree(base_image_folder)
        
os.mkdir(base_image_folder)
print("Wiped "+os.path.basename(base_image_folder))

#loop over all cases
for case_num in case_nums:
        print("case number "+str(case_num))
        
        case_num_folder = os.path.join(base_image_folder, "case_"+str(case_num))
        
        os.mkdir(case_num_folder)
        
        param_2 = 0
        if case_num == 0 or case_num == 4: param_2 = 1

        image_file_base = os.path.join(case_num_folder, "case_"+str(case_num)+"-"+str(param_2)+"_")


        #create list of parameters
        step_size = .01
        parameters_list = tuple(round(k,2) for k in np.arange(0, 1+step_size, step_size))

        A,B = tu.case(case_num)

        ####run in parallel on single node with multiprocessing####

        #function to run in parallel (makes this considerably faster)
        def generate_plot_for_k_value(k):
                imagefile = image_file_base+str(k)+".png"
                
                imagefile = image_file_base+str(k)+"-"+str(param_2)+".png"

                tu.plot_hausdorff_solution(A, param_2, B, 
                                    show_labels = True, 
                                    k=k,
                                    k_buffer = 0,
                                    imagefile = imagefile,
                                    dpi = 400)
                                    
                result = imageio.imread(imagefile)
                return result

        #get number of cpus available
        num_processes = mp.cpu_count()

        with mp.Pool(num_processes) as p:
                start = time.time()
                frames:list = p.map(generate_plot_for_k_value, parameters_list)
                
        ####run in parallel on single node with multiprocessing####

        #generate mp4 file
        MP4_PATH = image_file_base+"animation.mp4"
        fps=2
        writer = imageio.get_writer(MP4_PATH, fps=fps*5, macro_block_size=1)
        for frame in frames:
                writer.append_data(frame)
        writer.close()
        print("animation generated")
