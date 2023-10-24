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
out_main = '/global/cfs/cdirs/des/giannini/ggl/bfd_maglim_new2/'


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

# source redshift bins
# zs_bins = [
#            [0.20,0.43], 
#            [0.43,0.63], 
#            [0.63,0.90], 
#            [0.90,1.30],
#           ]
zs_bins = [
           [0.000, 0.405], 
           [0.405, 0.665], 
           [0.665, 0.960], 
           [0.96, 2.0],
          ]
# bin_edges =  [0.0, 0.405, 0.665, 0.96, 2.0]  
# {'bin0': 37687150, 'bin1': 37220369, 'bin2': 37766153, 'bin3': 37420019}


" Bins to run "
# the lens bins to run
l_bins = [0]
# l_bins = [0,1,2,3,4,5]
# the source bins to run
s_bins = [0,1,2,3]

source_cat = 'bfd'

" LSS weights "
# whether LSS weights will be applied to gamma_t
use_LSSweight = False

" Lens, source and random-point data files "
# lens galaxies
# data_lens = '/global/cfs/projectdirs/des/nweaverd/des_y6/maglim_plusplus/y6maglim_VIPERS_a18_b4_v2.0_pzcols_wisecols.fits'
data_lens = '/global/cfs/cdirs/des/y6kp-cats/maglim_2023-10-16.hdf5'
# data_lens = '/global/cfs/cdirs/des/giannini/mag_lim_lens_sample_combined_jointmask_sample_4ggl.fits'

# source galaxies
response = ['/global/cfs/cdirs/des/giannini/ggl/v5a_sompz/Response_bin0.txt',
            '/global/cfs/cdirs/des/giannini/ggl/v5a_sompz/Response_bin1.txt', 
            '/global/cfs/cdirs/des/giannini/ggl/v5a_sompz/Response_bin2.txt', 
            '/global/cfs/cdirs/des/giannini/ggl/v5a_sompz/Response_bin3.txt']

# source galaxies
# data_source = ['/global/cfs/cdirs/des/giannini/cats.yaml']
data_source = '/global/cfs/cdirs/des/y6kp-cats/metadetect_2023-10-16.hdf5'
data_source_bfd = ['/global/cfs/cdirs/des/giannini/ggl/v5a_sompz_jointmask/bfd_bin0.fits', 
                   '/global/cfs/cdirs/des/giannini/ggl/v5a_sompz_jointmask/bfd_bin1.fits',
                   '/global/cfs/cdirs/des/giannini/ggl/v5a_sompz_jointmask/bfd_bin2.fits',
                   '/global/cfs/cdirs/des/giannini/ggl/v5a_sompz_jointmask/bfd_bin3.fits']

# data_source = ['/global/cfs/cdirs/des/y6-shear-catalogs/Y6A2_METADETECT_V5a/metadetect_desdmv5a_cutsv5.h5']
# data_source = ['/global/cfs/cdirs/des/giannini/ggl/v5a_sompz/metadetect_bin0.fits',
               # '/global/cfs/cdirs/des/giannini/ggl/v5a_sompz/metadetect_bin1.fits', 
               # '/global/cfs/cdirs/des/giannini/ggl/v5a_sompz/metadetect_bin2.fits', 
               # '/global/cfs/cdirs/des/giannini/ggl/v5a_sompz/metadetect_bin3.fits']

# random points
# data_randoms = '/global/cfs/cdirs/des/giannini/Y3_maglim_randoms.fits'
data_randoms = '/global/cfs/cdirs/des/y6kp-cats/randoms_2023-10-16.hdf5'


" Boost factor "
use_boosts = True

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
bin_slop = 0.1 # supposed to be 0.0

# limits of theta in arcmins
theta_lims = [2.5, 250.]

# number of theta bins
ang_nbins = 20

# low memory - reduces memory usage, sacrificing speed
treecorr_low_mem = False


" Output directory for randoms and Jackknife patching "
path_out_rand = out_main+'/randoms_%dxlenses/'%(rand_fact)
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
