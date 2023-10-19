import numpy as np
import os
import sys
import yaml
import treecorr
import h5py as h5
from astropy.io import fits
import astropy.io.fits as pf
import healsparse
import joblib
import healpy as hp
from astropy.table import Table
import pandas as pd

path = os.path.abspath("./input/")
sys.path.append(path)

import params as par




def convert_to_pix_coord(ra, dec, nside=1024,nest=False):
    """
    Converts RA,DEC to hpix coordinates
    """

    theta = (90.0 - dec) * np.pi / 180.
    #print theta
    phi = ra * np.pi / 180.
    pix = hp.ang2pix(nside, theta, phi, nest=nest)

    return pix
        
# def prepare_source_metadetect(self, data_file, zs_bin=None):
#     """
#     Loads Y6 source galaxy data

#     """


    # met = h5.File(data_file)

print (par.zs_bins) 

met = h5.File('/global/cfs/cdirs/des/y6-shear-catalogs/Y6A2_METADETECT_V5a/metadetect_desdmv5a_cutsv5.h5')

ra_source = np.array(met['mdet/noshear']['ra'])
dec_source = np.array(met['mdet/noshear']['dec'])
pix_source = convert_to_pix_coord(ra_source,dec_source, nside=4096)

# load mask 1
hmap = healsparse.HealSparseMap.read('/global/cfs/cdirs/des/y6-shear-catalogs/y6-combined-hleda-gaiafull-des-stars-hsmap16384-nomdet-v3.fits')
# hmap = healsparse.HealSparseMap.read('/global/cfs/cdirs/des/y6-shear-catalogs/y6-combined-hleda-gaiafull-des-stars-hsmap16384-nomdet-v3.fits')
mask_1 = hmap.get_values_pos(ra_source, dec_source, valid_mask=True)

# load mask 2
lss_mask = pf.open('/global/cfs/projectdirs/des/monroy/Y6A2/maglim/sp_outliers_analysis/maglim_mask/jointmasks_fiducial/Y6LSSBAO_V2_MASK_WITHDEPTH_up_to_22.2_jointmask_3.5iqr_sps_0.01percent_sb_mean_0.5max_val_neb_mean_gcs_bit64_joint_vl05_vl10_vlim_zmax_gaia512_shear_flim0.8_bad_regions_gold20.fits.gz')
pix_lens_mask = lss_mask[1].data['HPIX_4096']
mask_2 = np.in1d(pix_source,pix_lens_mask)

joint_pix_mask_source = mask_1 & mask_2 

print ('done loading the masks')
print ('now i apply the masks')

sompz = h5.File('/global/cfs/cdirs/des/acampos/sompz_output/y6_data_preliminary/sompz_y6_data_preliminary.hdf5')
source = {}
source['ra']  = ra_source[joint_pix_mask_source]
source['dec'] = dec_source[joint_pix_mask_source]
source['e1']  = np.array(met['mdet/noshear']['gauss_g_1'])[joint_pix_mask_source]
print ('i am doing it!!')
source['e2']  = np.array(met['mdet/noshear']['gauss_g_2'])[joint_pix_mask_source]
source['g_cov_1_1']  = np.array(met['mdet/noshear']['gauss_g_cov_1_1'])[joint_pix_mask_source]
source['g_cov_2_2']  = np.array(met['mdet/noshear']['gauss_g_cov_2_2'])[joint_pix_mask_source]


source['som_bin'] = np.array(sompz['catalog/sompz/noshear/bhat'])[joint_pix_mask_source]

print ('now i look at the sheared cats')

ra_source_1p  = np.array(met['mdet/1p']['ra'])
dec_source_1p = np.array(met['mdet/1p']['dec'])
ra_source_2p  = np.array(met['mdet/2p']['ra'])
dec_source_2p = np.array(met['mdet/2p']['dec'])
ra_source_1m  = np.array(met['mdet/1m']['ra'])
dec_source_1m = np.array(met['mdet/1m']['dec'])
ra_source_2m  = np.array(met['mdet/2m']['ra'])
dec_source_2m = np.array(met['mdet/2m']['dec'])

shear_mask_shape_cat_1p = hmap.get_values_pos(ra_source_1p, dec_source_1p, valid_mask=True)
shear_mask_shape_cat_2p = hmap.get_values_pos(ra_source_2p, dec_source_2p, valid_mask=True)
shear_mask_shape_cat_1m = hmap.get_values_pos(ra_source_1m, dec_source_1m, valid_mask=True)
shear_mask_shape_cat_2m = hmap.get_values_pos(ra_source_2m, dec_source_2m, valid_mask=True)
nside = 4096 #16384
pix_source_1p = convert_to_pix_coord(ra_source_1p,dec_source_1p, nside=nside)
pix_source_2p = convert_to_pix_coord(ra_source_2p,dec_source_2p, nside=nside)
pix_source_1m = convert_to_pix_coord(ra_source_1m,dec_source_1m, nside=nside)
pix_source_2m = convert_to_pix_coord(ra_source_2m,dec_source_2m, nside=nside)

mask2_1p = np.in1d(pix_source_1p,pix_lens_mask)
mask2_2p = np.in1d(pix_source_2p,pix_lens_mask)
mask2_1m = np.in1d(pix_source_1m,pix_lens_mask)
mask2_2m = np.in1d(pix_source_2m,pix_lens_mask)


joint_pix_mask_source_1p = shear_mask_shape_cat_1p & mask2_1p
joint_pix_mask_source_2p = shear_mask_shape_cat_2p & mask2_2p
joint_pix_mask_source_1m = shear_mask_shape_cat_1m & mask2_1m
joint_pix_mask_source_2m = shear_mask_shape_cat_2m & mask2_2m


sombin1p = np.array(sompz['catalog/sompz/1p/bhat'])[joint_pix_mask_source_1p]
sombin2p = np.array(sompz['catalog/sompz/2p/bhat'])[joint_pix_mask_source_2p]
sombin1m = np.array(sompz['catalog/sompz/1m/bhat'])[joint_pix_mask_source_1m]
sombin2m = np.array(sompz['catalog/sompz/2m/bhat'])[joint_pix_mask_source_2m]


print (par.zs_bins) 

for zs_bin in range(len(par.zs_bins)):
    
    print ('zs_bin', zs_bin)
    
    maskz = (source['som_bin'] == zs_bin) 

    ra_s = source['ra'][maskz]
    dec_s = source['dec'][maskz]
    e1_s = source['e1'][maskz]
    e2_s = source['e2'][maskz]
    g_cov_1_1 = source['g_cov_1_1'][maskz]
    g_cov_2_2 = source['g_cov_2_2'][maskz]

    maskz1p = (sombin1p == zs_bin) 
    maskz2p = (sombin2p == zs_bin) 
    maskz1m = (sombin1m == zs_bin) 
    maskz2m = (sombin2m == zs_bin) 

    
    g1_1p = np.array(met['mdet/1p/gauss_g_1'])
    g1_2p = np.array(met['mdet/2p/gauss_g_1'])
    g1_1m = np.array(met['mdet/1m/gauss_g_1'])
    g1_2m = np.array(met['mdet/2m/gauss_g_1'])
    g2_1p = np.array(met['mdet/1p/gauss_g_2'])
    g2_2p = np.array(met['mdet/2p/gauss_g_2'])
    g2_1m = np.array(met['mdet/1m/gauss_g_2'])
    g2_2m = np.array(met['mdet/2m/gauss_g_2'])
    g1_noshear = np.array(met['mdet/noshear/gauss_g_1'])
    g2_noshear = np.array(met['mdet/noshear/gauss_g_2'])
  

    def get_shear_weights(cols, weights_file, weight_type):

        def _assign_loggrid(x, y, xmin, xmax, xsteps, ymin, ymax, ysteps):
            from math import log10
            # return x and y indices of data (x,y) on a log-spaced grid that runs from [xy]min to [xy]max in [xy]steps

            logstepx = log10(xmax/xmin)/xsteps
            logstepy = log10(ymax/ymin)/ysteps

            indexx = (np.log10(x/xmin)/logstepx).astype(int)
            indexy = (np.log10(y/ymin)/logstepy).astype(int)

            indexx = np.maximum(indexx,0)
            indexx = np.minimum(indexx, xsteps-1)
            indexy = np.maximum(indexy,0)
            indexy = np.minimum(indexy, ysteps-1)

            return indexx,indexy

        def _find_shear_weight(d, wgt_dict, snmin, snmax, sizemin, sizemax, steps):

            if wgt_dict is None:
                weights = np.ones(len(d))
                return weights

            shear_wgt = wgt_dict['weight']
            indexx, indexy = _assign_loggrid(d['s2n'], d['T_ratio'], snmin, snmax, steps, sizemin, sizemax, steps)
            weights = np.array([shear_wgt[x, y] for x, y in zip(indexx, indexy)])

            return weights

        #import pdb ; pdb.set_trace()
        if weight_type == 's2n_sizer':
            # pickle file that defines w(S/N, size)
            with open(weights_file, 'rb') as handle:
                wgt_dict = pickle.load(handle)
            ## TO-DO: make snmin, snmax, sizemin, sizemax available in config file. 
            shear_wgt = _find_shear_weight(cols, wgt_dict, 10, 1000, 0.5, 5.0, 20)
        elif weight_type == 'shape_err':
            shear_wgt = 1/(0.17**2 + 0.5*(cols['g_cov_1_1'] + cols['g_cov_2_2']))

        shear_wgt[np.isnan(shear_wgt)] = 0.

        return shear_wgt

    weights_file = '/global/cfs/cdirs/des/myamamot/y6_shear_catalogs/Y6A2_METADETECT_V4/inverse_variance_weight_v3_s2n_10-300_Tratio_0.5-5.pickle'
    weight_type = 'shape_err'

    source['w'] = get_shear_weights(source, weights_file, weight_type)
    w_g = source['w'][maskz]
    print ('wg', w_g[:10])
    
      
    mean_g1 = np.sum(g1_noshear[joint_pix_mask_source][maskz])/len(g1_noshear[joint_pix_mask_source][maskz])
    mean_g1p = np.sum(g1_1p[joint_pix_mask_source_1p][maskz1p])/len(g1_1p[joint_pix_mask_source_1p][maskz1p])
    mean_g1m = np.sum(g1_1m[joint_pix_mask_source_1m][maskz1m])/len(g1_1m[joint_pix_mask_source_1m][maskz1m])
    R11 = (mean_g1p-mean_g1m)/2/0.01

        
    mean_g2 = np.sum(g2_noshear[joint_pix_mask_source][maskz])/len(g2_noshear[joint_pix_mask_source][maskz])
    mean_g2p = np.sum(g2_2p[joint_pix_mask_source_2p][maskz2p])/len(g2_2p[joint_pix_mask_source_2p][maskz2p])
    mean_g2m = np.sum(g2_2m[joint_pix_mask_source_2m][maskz2m])/len(g2_2m[joint_pix_mask_source_2m][maskz2m])
    R22 = (mean_g2p-mean_g2m)/2/0.01
    R_g = (R11 + R22)/2
    
    
#     mean_g1 = np.sum(g1_noshear[joint_pix_mask_source][maskz])/len(g1_noshear[joint_pix_mask_source][maskz])
#     mean_g1p = np.sum(g1_1p[joint_pix_mask_source_1p][maskz1p])/len(g1_1p[joint_pix_mask_source_1p][maskz1p])
#     mean_g1m = np.sum(g1_1m[joint_pix_mask_source_1m][maskz1m])/len(g1_1m[joint_pix_mask_source_1m][maskz1m])
#     R11 = (mean_g1p-mean_g1m)/2/0.01

        
#     mean_g2 = np.sum(g2_noshear[joint_pix_mask_source][maskz])/len(g2_noshear[joint_pix_mask_source][maskz])
#     mean_g2p = np.sum(g2_2p[joint_pix_mask_source_2p][maskz2p])/len(g2_2p[joint_pix_mask_source_2p][maskz2p])
#     mean_g2m = np.sum(g2_2m[joint_pix_mask_source_2m][maskz2m])/len(g2_2m[joint_pix_mask_source_2m][maskz2m])
#     R22 = (mean_g2p-mean_g2m)/2/0.01
#     R_g = (R11 + R22)/2


    mm = pd.DataFrame()

    mm['ra_s'] = ra_s 
    mm['dec_s'] = dec_s
    mm['e1_s'] = e1_s 
    mm['e2_s'] = e2_s 
    mm['g_cov_1_1'] = g_cov_1_1
    mm['g_cov_2_2'] = g_cov_2_2
    mm['w_g'] = w_g
    # mm['R_g'] = R_g


    # path = '/global/cfs/cdirs/des/giannini/ggl/metad_maglim_fidmask/'
    path = par.out_main
    if not os.path.exists(path):
        os.makedirs(path)

    meta_table = Table.from_pandas(mm)
    meta_table.write(path+'/metadetect_bin{}.fits'.format(zs_bin), overwrite = True)

    response = np.array([R_g, R11, R22])
    np.savetxt(path+'/Response_bin{}.txt'.format(zs_bin), response)


    # return source['ra'] , source['dec'], source['e1'], source['e2'], R_g, w_g
    # return ra_s, dec_s, e1_s, e2_s, R_g, w_g



