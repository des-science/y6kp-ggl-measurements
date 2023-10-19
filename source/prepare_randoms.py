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
import joblib
import healsparse


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


print (par.zs_bins) 




# load into memory randoms
randoms = pf.open('/global/cfs/cdirs/des/jelvinpo/y6_randoms/y6maglim_randoms_rr40.0_y6_maglim_jointmask_fiducial.fits')
ra_randoms = np.array(randoms[1].data['RA'])
dec_randoms = np.array(randoms[1].data['DEC'])
nside = 4096
pix_randoms = convert_to_pix_coord(ra_randoms,dec_randoms, nside=nside)

# load mask 1:
hmap = healsparse.HealSparseMap.read('/global/cfs/cdirs/des/y6-shear-catalogs/y6-combined-hleda-gaiafull-des-stars-hsmap16384-nomdet-v3.fits')
# y6-combined-hleda-gaiafull-des-stars-hpmap4096-nomdet-fracdet-v3.fits

mask_1 = hmap.get_values_pos(ra_randoms,dec_randoms, valid_mask=True)

# load mask_2
lss_mask = pf.open('/global/cfs/projectdirs/des/monroy/Y6A2/maglim/sp_outliers_analysis/maglim_mask/FIDUCIAL_MASK/y6_maglim_jointmask_fiducial.fits.gz')
# lss_mask = pf.open('/global/cfs/projectdirs/des/monroy/Y6A2/maglim/sp_outliers_analysis/maglim_mask/jointmasks_fiducial/Y6LSSBAO_V2_MASK_WITHDEPTH_up_to_22.2_jointmask_3.5iqr_sps_0.01percent_sb_mean_0.5max_val_neb_mean_gcs_bit64_joint_vl05_vl10_vlim_zmax_gaia512_shear_flim0.8_bad_regions_gold20.fits.gz')
pix_lens_mask = lss_mask[1].data['HPIX_4096']
mask_2 = np.in1d(pix_randoms,pix_lens_mask)

joint_pix_mask_randoms = mask_1 & mask_2 

print ('done loading the masks')
print ('now i apply the masks')

from astropy.table import Table
import pandas as pd


mm = pd.DataFrame()

mm['RA'] = ra_randoms[joint_pix_mask_randoms]
mm['DEC'] = dec_randoms[joint_pix_mask_randoms]


path = par.out_main
if not os.path.exists(path):
    os.makedirs(path)

stupid = Table.from_pandas(mm)
stupid.write(path+'/randoms++.fits', overwrite = True)
        
    
    