"""
----------------------------
Main input file for Y6KP GGL
----------------------------
"""
import os
import sys

here = os.path.abspath(__file__)
sys.path.append(here)

" Define redshift bins for lens and source galaxies "
# lens redsift bins
zl_bins = [
           [0.20,0.40], 
           [0.40,0.55], 
           [0.55,0.70], 
           [0.70,0.85],
          ]

# source redshift bins
zs_bins = [
           [0.20,0.43], 
           [0.43,0.63], 
           [0.63,0.90], 
           [0.90,1.30],
          ]

" Bins to run "
l_bins  = [0,1,2,3]
s_bins  = [0,1,2,3]

" Mask file of lenses "
# dictionary with masks for all lens bins;
# must be 0's and 1's in ".npy" or ".npz" format;
# 'key' convention in dictionary: 'l{lbin}m{mbin}'
lens_mask_file = '...path_to_mask...'
lens_mask_nested = False

" Mask file of randoms "
# dictionary with masks for all lens bins;
# must be 0's and 1's in ".npy" or ".npz" format;
# 'key' convention in dictionary: 'l{lbin}m{mbin}'
randoms_mask_file = '...path_to_mask...'
randoms_mask_nested = False

" LSS weights "
# mask need to be a dictionary, 
# with keys "weightLSS_name_zbin_mbin", where zbin and mbins are integers,
# in ".npz" format
use_LSSweight = False
weightLSS_file = '...path_to_file_with_list_of_weights...'
weightLSS_nest = False

" Current directory "
here = os.path.dirname(os.path.abspath(__file__))

" Jackknife's source code directory "
source_dir = './source/'

" Lens and source data files "
# lens galaxies; structure: [...,[lens file for bin i],...]
lens_dir         = '...path_to_folder_with_lens_data_files...'

# source galaxies
# note: 'data_source_galaxies' can be one of ['metacal_5sels', 'mice']
dir_descat  = '/project2/chihway/data/des_y3_catalogs/y3kp_sample/Feb2021/'
file_yaml   = 'cats.yaml'
data_source = [dir_descat+file_yaml,
               dir_descat+file_yaml,
               dir_descat+file_yaml,
               dir_descat+file_yaml]
data_source_galaxies = 'metacal_5sels'

" Boost factor "
use_boosts  = True

" Random point subtraction "
# use random subtraction
use_randoms = True
# factor of randoms/lenses
rand_fact = 30

" Treecorr settings "
# number of Jackknife regions
n_jck = 150

# resolution of grid (make it a power of 2)
nside = 4096

# allowed bin slop for treecorr
bin_slop = 0.0

# limits of theta in arcmins
theta_lims = [0.25, 250.]

# number of theta bins
ang_nbins = 30

# look for nearby sources to run treecorr
# for each of the Jackknife patches
look_for_neighbors = True

# resolution to find sources nearby lens patch
nside_nearby = 4

# low memory - reduces memory usage, sacrificing speed
treecorr_low_mem = False

" Method for Jackknife calculation "
# available methods:
#   I)  'remove': remove the area and keep rest of sky
#   II) 'keep'  : keep just that one area and remove rest of sky
jck_method = 'keep'

" Main output folder "
out_main = './Output/'

" Output directory for randoms and Jackknife patching "
path_out_rand = out_main+'/randoms_%dxlenses/'%(rand_fact)
path_out_JK   = out_main+'/Jackknife/'
path_JK_info  = path_out_JK+'/info/'

" Output folder for boost factors "
path_out_boost = out_main+'/boosts/'
path_JK_boost  = path_out_JK+'/boosts_JK/'
path_JK_cov_bf = path_out_JK+'/boosts_covariance/'

" Output directories for gamma_t "
out_main_gt       = out_main
path_out_gt       = out_main_gt+'/gammat/'
path_out_gt_rand  = out_main_gt+'/gammat_randoms/'
path_out_extra_gt = out_main_gt+'/gammat_extra/'
path_out_gx       = out_main_gt+'/gammax/'
path_out_gx_rand  = out_main_gt+'/gammax_randoms/'
path_JK_cov_gx    = path_out_JK+'/gammax_covariance/'
#
path_JK_gt       = path_out_JK+'/gammat_JK/'
path_JK_gx       = path_out_JK+'/gammax_JK/'
path_JK_rand     = path_out_JK+'/gammat_randoms_JK/'
path_JK_extra_gt = path_out_JK+'/gammat_extra_JK/'
path_JK_cov_gt   = path_out_JK+'/gammat_covariance/'
path_out_shot_gt = out_main_gt+'/gammat_shot_noise/'

" Output directories for wtheta "
out_main_wt       = out_main
path_out_wt       = out_main_wt+'/wtheta/'
path_out_extra_wt = out_main_wt+'/wtheta_extra/'
#
path_JK_cov_wt   = path_out_JK+'/wtheta_covariance/'
path_JK_wt       = path_out_JK+'/wtheta_JK/'
path_JK_extra_wt = path_out_JK+'/wtheta_extra_JK/'
path_out_shot_wt = out_main_wt+'/wtheta_shot_noise/'
