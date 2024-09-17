"""
----------------------------
Main input file for Y6KP GGL
----------------------------
"""
import os
import sys

here = os.path.abspath(__file__)
sys.path.append(here)

" Main output folder "
out_main = '/global/cfs/projectdirs/des/awhyley/ggl/fullscale_070624'

" Bins to run "
# the lens bins to run
l_bins = [0,1,2,3,4,5]
#l_bins = [0]

# the source bins to run
s_bins = [0,1,2,3]
# s_bins = [0]

" Source catalog "
# bdf or metadetect
source_cat = 'metadetect' 


" Define redshift bins for lens and source galaxies "
# lens redsift bins
zl_bins = [
            [0.20, 0.40], 
            [0.40, 0.55], 
            [0.55, 0.70], 
            [0.70, 0.85],
            [0.85, 0.95],
            [0.95, 1.05],
          ]
#zl_bins = [[0.20, 1.05]]

# source redshift bins
zs_bins = [
           [0.000, 0.405], 
           [0.405, 0.665], 
           [0.665, 0.960], 
           [0.96, 2.0],
          ]

" Lens, source and random-point data files "
# lens galaxies
data_lens = '/global/cfs/cdirs/des/y6kp-cats/2024-07-15/desy6kp_cats_2024-07-15.hdf5'

# random points
data_randoms = '/global/cfs/cdirs/des/y6kp-cats/2024-07-15/desy6kp_cats_2024-07-15.hdf5'


# source galaxies
data_source = '/global/cfs/cdirs/des/giannini/ggl/data/2024-08-26/metadetect_v6_UNBLINDED_2024-08-26.hdf5'
response = ['/global/cfs/cdirs/des/giannini/ggl/mdet_2024-08-26_response_bin1.txt',
            '/global/cfs/cdirs/des/giannini/ggl/mdet_2024-08-26_response_bin2.txt', 
            '/global/cfs/cdirs/des/giannini/ggl/mdet_2024-08-26_response_bin3.txt', 
            '/global/cfs/cdirs/des/giannini/ggl/mdet_2024-08-26_response_bin4.txt']

#Change for each band and size/shape for PSF test

#this is old but without it it does not run
data_source_bfd = '/global/cfs/cdirs/des/y6kp-cats/2023-10-16/desy6kp_cats_2023-10-16.hdf5'


# Input datavector for substitution of ggl measurement
dv_input = '/global/cfs/cdirs/des/giannini/blinding/data_vectors/2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate_nosourcesreals_sompzmean.fits'

# Output datavector 
dv_output = '/global/cfs/projectdirs/des/awhyley/ggl/fullscale_070624/fullscale_dv.fits' 
# dv_output = '/global/cfs/cdirs/des/giannini/y6-3x2pt/blinding/data_vectors/2pt_y3dv_y6_ggl_maglim_metadet__2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate_nosourcesreals_sompzmean_binslop0.fits'
#Doesn't matter for psf test



" LSS weights "
# whether LSS weights will be applied to gamma_t
use_LSSweight = True

" Shear weights "
# whether shear weights will be applied to gamma_t
use_shearweight = False

" Boost factor "
use_boosts = True

" Reponse "
use_response = True

" Random point subtraction "
# use random subtraction
use_randoms = True
# factor of randoms/lenses
# rand_fact = 30

" Treecorr settings "
# number of Jackknife regions
n_jck = 250

# resolution of grid (make it a power of 2)
nside = 4096

# allowed bin slop for treecorr
bin_slop = 0.1 # supposed to be 0.0
angle_slop = 0.05

# limits of theta in arcmins
theta_lims = [2.5, 995.26792638]

# number of theta bins
# ang_nbins = 20
ang_nbins = 26

# low memory - reduces memory usage, sacrificing speed
treecorr_low_mem = False

# used for calculating scale dependant response
calc_scale_dependant_response = False


" Output directory for randoms and Jackknife patching "
path_out_rand = out_main+'/randoms_lenses/'
path_out_JK = out_main+'/Jackknife/'
path_JK_info = path_out_JK+'/info/'

" Output folder for boost factors "
path_out_boost = out_main+'/boosts/'
path_JK_boost = path_out_JK+'/boosts_JK/'
path_JK_cov_bf = path_out_JK+'/boosts_covariance/'

" Output directories for gamma_t "
out_main_gt = out_main
path_out_gt = out_main_gt+'/gammat/'
path_out_gt_rand = out_main_gt+'/gammat_randoms/'
path_out_extra_gt = out_main_gt+'/gammat_extra/'
path_out_gx = out_main_gt+'/gammax/'
path_out_gx_rand = out_main_gt+'/gammax_randoms/'
path_JK_cov_gx = path_out_JK+'/gammax_covariance/'
#
path_JK_gt = path_out_JK+'/gammat_JK/'
path_JK_gx = path_out_JK+'/gammax_JK/'
path_JK_rand = path_out_JK+'/gammat_randoms_JK/'
path_JK_extra_gt = path_out_JK+'/gammat_extra_JK/'
path_JK_cov_gt = path_out_JK+'/gammat_covariance/'
path_JK_cov_gt_rand = path_out_JK+'/gammat_randoms_covariance/'
path_out_shot_gt = out_main_gt+'/gammat_shot_noise/'
