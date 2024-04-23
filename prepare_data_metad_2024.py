import numpy as np
import pickle
import os
import shutil
import pandas as pd
import sys
import yaml
import astropy.io.fits as pf
from astropy.table import Table
import healpy as hp
import h5py as h5
import healsparse
import joblib
import healpy as hp

def convert_to_pix_coord(ra, dec, nside=1024,nest=False):
    """
    Converts RA,DEC to hpix coordinates
    """

    theta = (90.0 - dec) * np.pi / 180.
    #print theta
    phi = ra * np.pi / 180.
    pix = hp.ang2pix(nside, theta, phi, nest=nest)

    return pix

def IndexToDeclRa(index, nside=1024,nest= False):
    theta,phi=hp.pixelfunc.pix2ang(nside ,index,nest=nest)
    return -np.degrees(theta-np.pi/2.),np.degrees(phi)

def generate_randoms_radec(minra, maxra, mindec, maxdec, Ngen, raoffset=0):
    r = 1.0
    # this z is not redshift!
    zmin = r * np.sin(np.pi * mindec / 180.)
    zmax = r * np.sin(np.pi * maxdec / 180.)
   
    # parity transform from usual, but let's not worry about that
    phimin = np.pi / 180. * (minra - 180 + raoffset)
    phimax = np.pi / 180. * (maxra - 180 + raoffset)
   
    # generate ra and dec
    z_coord = np.random.uniform(zmin, zmax, Ngen)  # not redshift!
    phi = np.random.uniform(phimin, phimax, Ngen)
    dec_rad = np.arcsin(z_coord / r)
   
    # convert to ra and dec
    ra = phi * 180 / np.pi + 180 - raoffset
    dec = dec_rad * 180 / np.pi
    return ra, dec


y6kp = h5.File('/global/cfs/cdirs/des/y6kp-cats/2023-10-16/desy6kp_cats_2023-10-16.hdf5')

y6kp['desy6kp/mdet/noshear/

# load catalog
met = h5.File('/global/cfs/cdirs/des/y6-shear-catalogs/Y6A2_METADETECT_V5a/metadetect_desdmv5a_cutsv5.h5')

ra_source = np.array(met['mdet/noshear']['ra'])
dec_source = np.array(met['mdet/noshear']['dec'])
nside = 4096
pix_source = convert_to_pix_coord(ra_source,dec_source, nside=nside)


# load mask 1
hmap = healsparse.HealSparseMap.read('/global/cfs/cdirs/des/y6-shear-catalogs/y6-combined-hleda-gaiafull-des-stars-hsmap16384-nomdet-v3.fits')
mask_1 = hmap.get_values_pos(ra_source, dec_source, valid_mask=True)

# load mask 2
lss_mask = pf.open('/global/cfs/projectdirs/des/monroy/Y6A2/maglim/sp_outliers_analysis/maglim_mask/FIDUCIAL_MASK/y6_maglim_jointmask_fiducial.fits.gz')
# lss_mask = pf.open('/global/cfs/projectdirs/des/monroy/Y6A2/maglim/sp_outliers_analysis/maglim_mask/jointmasks_fiducial/Y6LSSBAO_V2_MASK_WITHDEPTH_up_to_22.2_jointmask_3.5iqr_sps_0.01percent_sb_mean_0.5max_val_neb_mean_gcs_bit64_joint_vl05_vl10_vlim_zmax_gaia512_shear_flim0.8_bad_regions_gold20.fits.gz')
pix_lens_mask = lss_mask[1].data['HPIX_4096']
mask_2 = np.in1d(pix_source,pix_lens_mask)

joint_pix_mask_source = mask_1 & mask_2 


sompz = h5.File('/global/cfs/cdirs/des/acampos/sompz_output/y6_data_preliminary/sompz_y6_data_preliminary.hdf5')
source = {}
source['ra']  = np.array(met['mdet/noshear']['ra'][joint_pix_mask_source])
source['dec'] = np.array(met['mdet/noshear']['dec'][joint_pix_mask_source])
source['e1']  = np.array(met['mdet/noshear']['gauss_g_1'][joint_pix_mask_source])
source['e2']  = np.array(met['mdet/noshear']['gauss_g_2'][joint_pix_mask_source])
source['g_cov_1_1']  = np.array(met['mdet/noshear']['gauss_g_cov_1_1'][joint_pix_mask_source])
source['g_cov_2_2']  = np.array(met['mdet/noshear']['gauss_g_cov_2_2'][joint_pix_mask_source])


source['som_bin'] = np.array(sompz['catalog/sompz/noshear/bhat'][joint_pix_mask_source])
print ('------> 0')


ra_source_1p  = np.array(met['mdet/1p']['ra'])
dec_source_1p = np.array(met['mdet/1p']['dec'])
ra_source_2p  = np.array(met['mdet/2p']['ra'])
dec_source_2p = np.array(met['mdet/2p']['dec'])
ra_source_1m  = np.array(met['mdet/1m']['ra'])
dec_source_1m = np.array(met['mdet/1m']['dec'])
ra_source_2m  = np.array(met['mdet/2m']['ra'])
dec_source_2m = np.array(met['mdet/2m']['dec'])


print ('------> 1')

shear_mask_shape_cat_1p = hmap.get_values_pos(ra_source_1p, dec_source_1p, valid_mask=True)
shear_mask_shape_cat_2p = hmap.get_values_pos(ra_source_2p, dec_source_2p, valid_mask=True)
shear_mask_shape_cat_1m = hmap.get_values_pos(ra_source_1m, dec_source_1m, valid_mask=True)
shear_mask_shape_cat_2m = hmap.get_values_pos(ra_source_2m, dec_source_2m, valid_mask=True)
nside = 4096 #16384
pix_source_1p = convert_to_pix_coord(ra_source_1p,dec_source_1p, nside=nside)
pix_source_2p = convert_to_pix_coord(ra_source_2p,dec_source_2p, nside=nside)
pix_source_1m = convert_to_pix_coord(ra_source_1m,dec_source_1m, nside=nside)
pix_source_2m = convert_to_pix_coord(ra_source_2m,dec_source_2m, nside=nside)

print ('------> 2')

mask2_1p = np.in1d(pix_source_1p,pix_lens_mask)
mask2_2p = np.in1d(pix_source_2p,pix_lens_mask)
mask2_1m = np.in1d(pix_source_1m,pix_lens_mask)
mask2_2m = np.in1d(pix_source_2m,pix_lens_mask)


joint_pix_mask_source_1p = shear_mask_shape_cat_1p & mask2_1p
joint_pix_mask_source_2p = shear_mask_shape_cat_2p & mask2_2p
joint_pix_mask_source_1m = shear_mask_shape_cat_1m & mask2_1m
joint_pix_mask_source_2m = shear_mask_shape_cat_2m & mask2_2m

print ('------> 3')

sombin1p = np.array(sompz['catalog/sompz/1p/bhat'][joint_pix_mask_source_1p])
sombin2p = np.array(sompz['catalog/sompz/2p/bhat'][joint_pix_mask_source_2p])
sombin1m = np.array(sompz['catalog/sompz/1m/bhat'][joint_pix_mask_source_1m])
sombin2m = np.array(sompz['catalog/sompz/2m/bhat'][joint_pix_mask_source_2m])




for zs_bin in range(4):

    print (zs_bin)
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

    mean_g1 = np.sum(met['mdet/noshear/gauss_g_1'][joint_pix_mask_source][maskz])/len(met['mdet/noshear/gauss_g_1'][joint_pix_mask_source][maskz])
    mean_g1p = np.sum(met['mdet/1p/gauss_g_1'][joint_pix_mask_source_1p][maskz1p])/len(met['mdet/1p/gauss_g_1'][joint_pix_mask_source_1p][maskz1p])
    mean_g1m = np.sum(met['mdet/1m/gauss_g_1'][maskz1m])/len(met['mdet/1m/gauss_g_1'][maskz1m])
    R11 = (mean_g1p-mean_g1m)/2/0.01

    mean_g2 = np.sum(met['mdet/noshear/gauss_g_2'][maskz])/len(met['mdet/noshear/gauss_g_2'][maskz])
    mean_g2p = np.sum(met['mdet/2p/gauss_g_2'][joint_pix_mask_source_2p][maskz2p])/len(met['mdet/2p/gauss_g_2'][joint_pix_mask_source_2p][maskz2p])
    mean_g2m = np.sum(met['mdet/2m/gauss_g_2'][joint_pix_mask_source_2m][maskz2m])/len(met['mdet/2m/gauss_g_2'][joint_pix_mask_source_2m][maskz2m])
    R22 = (mean_g2p-mean_g2m)/2/0.01
    R_g = (R11 + R22)/2
    
    
    from astropy.table import Table
    import pandas as pd


    mm = pd.DataFrame()
    
    mm['ra_s'] = ra_s 
    mm['dec_s'] = dec_s
    mm['e1_s'] = e1_s 
    mm['e2_s'] = e2_s 
    mm['g_cov_1_1'] = g_cov_1_1
    mm['g_cov_2_2'] = g_cov_2_2


    path = '/global/cfs/cdirs/des/giannini/ggl/metad_maglim_fidmask/'
    if not os.path.exists(path):
        os.makedirs(path)
    
    stupid = Table.from_pandas(mm)
    stupid.write(path+'/metadetect_bin{}.fits'.format(zs_bin), overwrite = True)
    
    response = np.array([R_g, R11, R22])
    np.savetxt(path+'/Response_bin{}.txt'.format(zs_bin), response)
    
    
    