"""
------------------------
Run the GGL measurements
------------------------
"""
" Preliminary setup "
import sys
import os

# parameters file
import info as info

# folder containing important stuff
sys.path.append(info.input_dir)

# code to run measurements
import ggly6 as ggl_run
run = ggl_run.GGL(input_dir=info.input_dir)
print ('info.input_dir', info.input_dir)

#------------------------------------------------#

"""
*******************
Do the calculations
*******************
"""
if __name__ == "__main__":

    # run code to generate measurements and covariance
    if info.NG_mode=='gt_only':
        print('Running gamma_t calculations (no Jackknife covariance)...')
        run.run_gammat()
        
    elif info.NG_mode=='gt_and_cov_parallel':
        print('Running gamma_t and Jackknife covariance calculations with MPI...')
        run.run_gammat_and_cov_parallel()
        
    elif info.NG_mode=='gt_and_cov':
        print('Running gamma_t and Jackknife covariance calculations...')
        run.run_gammat_and_cov()
        
    elif info.NG_mode=='NK':
        print('Running NK correlations for scale-dependant response test...')
        run.run_nk()

    elif info.NG_mode=='gt_and_cov_catalogsonly':
        print('Saving treecorr catalogs only...')
        run.run_gammat_and_cov_catalogsonly()        

        
    else:
        errmsg = "!!!Error: Requested NG mode not implemented"
        raise Exception(errmsg)

    exit()
