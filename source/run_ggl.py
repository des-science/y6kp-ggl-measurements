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

# code to run measurements
import ggly6 as ggl_run
run = ggl_run.GGL(input_dir=info.input_dir)

#------------------------------------------------#

"""
*******************
Do the calculations
*******************
"""
if __name__ == "__main__":

    # run code to generate measurements and covariance
    for corr in par.correlations:
        if corr=='NG':
            if par.NG_setting['corr_only']:
                print('Running gamma_t calculations...')
                jack_run.run_gammat()
            if par.NG_setting['corr_and_cov']:
                print('Running gamma_t and covariance calculations using TreeCorr only...')
                jack_run.run_gammat_and_cov()
            if par.NG_setting['corr_and_cov_myJK']:
                print('Running gamma_t and covariance calculations using my own code...')
                jack_run.run_jck_gammat()
            if par.NG_setting['postprocess']:
                print('Postprocessing outputs to calculate the derived quantities and covariances...')
                jack_run.postprocess_jackknife()
        elif corr=='NN':
            if par.NN_setting['corr_only']:
                print('Running wtheta calculations...')
                jack_run.run_wtheta()
            if par.NN_setting['corr_and_cov']:
                print('Running wtheta and covariance calculations using TreeCorr only...')
                jack_run.run_wtheta_and_cov()
        else:
            errmsg = '!!!Error: Correlation "{}" is not implemented, choose one of {}'.format(corr,par.correlations)
            raise Exception(errmsg)

    exit()
