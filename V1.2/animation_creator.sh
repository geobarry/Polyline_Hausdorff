#!/bin/bash

animation_creator()
{(

        set -Eeo pipefail
        
        local outputs_folder=$1
        local cwd=$2
        
        #set directory, check if it exists
        [ -z $cwd ] && cwd=$PWD
        [ ! -d $cwd ] && echo $cwd" does not exist" && exit 1
	
        cd $cwd
	
	#switch to environment
        eval "$(conda shell.bash hook)"
        conda activate s3plot_camerer
        
        #run script
        python -W ignore animation_creator.py
)}

animation_creator '/home/acamerer/Desktop/My_Project/Workflows/hausdroff_animation/outputs'
