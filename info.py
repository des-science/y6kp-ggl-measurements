"""
File containing the input directory and mode 
to start running the Y6 GGL measurements pipeline
"""
# path the the folder containing all the input files
# input_dir = '...path_to_input_dir...'
input_dir = './2025-02-03_fid_matchweight'

# mode to run NG correlations
# options:
# I)  'gt_only': only run the correlation funtions without covariance
# II) 'gt_and_cov': run gamma_t measurements with its Jackknife covariance
# III) 'gt_and_cov_parallel': run gamma_t measurements with its Jackknife covariance, using MPI -> srun --nodes=4 --tasks-per-node=1 python run_ggl.py
# IV) 'NK': run NK correlations for scale-dependant response test
NG_mode = 'gt_and_cov_parallel'
