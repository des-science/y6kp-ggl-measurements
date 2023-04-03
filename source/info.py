"""
File containing the necessary information to setup and 
start running the Y6 GGL measurements pipeline
"""
import os
import sys

here = os.path.abspath(__file__)
sys.path.append(here)

# path the the folder containing all the input files
input_dir = '...path_to_input_dir...'

# mode to run NG correlations
# options:
# I)  'gt_only': only run the correlation funtions without covariance
# II) 'gt_and_cov': run gamma_t measurements with its Jackknife covariance
NG_mode = 'gt_and_cov'
