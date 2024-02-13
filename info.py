"""
File containing the input directory and mode 
to start running the Y6 GGL measurements pipeline
"""
# path the the folder containing all the input files
# input_dir = '...path_to_input_dir...'
input_dir = './input'

# mode to run NG correlations
# options:
# I)  'gt_only': only run the correlation funtions without covariance
# II) 'gt_and_cov': run gamma_t measurements with its Jackknife covariance
# III) 'gt_and_cov_par': run gamma_t measurements with its Jackknife covariance, using MPI -> run using srun -n 24 python run_ggl.py
NG_mode = 'gt_and_cov'
