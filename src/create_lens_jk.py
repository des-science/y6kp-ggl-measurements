import numpy as np
import sys
import yaml
sys.path.append('../../destest/')
sys.path.append('../../kmeans_radec/')
import destest
#import pyfits as pf
import astropy.io.fits as fits
from info import zbins, config, paths
import destest_functions
import kmeans_radec
import os
import jk
import ipdb

print config

def create_lens(lens_selector, ran_selector, pz_selector, gold_selector, zbins, path_save, sample):

        """
        Function that saves lens.fits and random.fits files, after splitting into JK regions.
        lens_selector: destest object to load the lens sample from the hdf5 file.
        ran_selector: destest object to load the random sample from the hdf5 file.
        pz_selector: destest object to load the redshift sample from the hdf5 file, 
                     only necessary for the alternative lens sample. Redmagic has their own z already.
        gold_selector: only necessary for maglimit sample.
        zbins: dictionary from info with the limits of the sample. Here we only use the low and high limits 
               to avoid selecting galaxies we will not use later.
        path_save: where the files are saved.
        sample: string, can be 'redmagic' or 'maglim'.
        """
        print 'Sample:', sample
        print 'Saving lens sample here:', path_save

        ra_l  = lens_selector.get_col('ra')[0]
        dec_l = lens_selector.get_col('dec')[0]
        if sample == 'redmagic': 
                z_l = lens_selector.get_col('zredmagic')[0]
                zerr_l = lens_selector.get_col('zredmagic_e')[0]
                ids = lens_selector.get_col('coadd_object_id')[0]
        if sample == 'maglim':
                z_l = pz_selector.get_col('zmean_sof')[0]
                zerr_l = pz_selector.get_col('zmc_sof')[0]
                ids = gold_selector.get_col('coadd_object_id')[0]
        w_l = lens_selector.get_col('weight')[0]
        print 'Weights lenses:', w_l
        assert len(ra_l)==len(ids), 'Something is wrong.'
        print 'Range of lens weights:', np.min(w_l), np.max(w_l)

        # Select galaxies that belong to some lens z-bin
        zbins = zbins['lims']
        print 'zbins:', zbins
        maskzl = ((z_l>zbins[0])&(z_l<zbins[-1]))
        ra_l = ra_l[maskzl]
        dec_l = dec_l[maskzl]
        z_l = z_l[maskzl]
        zerr_l = zerr_l[maskzl]
        ids = ids[maskzl]
        w_l = w_l[maskzl]

        #Load ra,dec, from random catalog 
        ra_r = ran_selector.get_col('ra')[0] 
        dec_r = ran_selector.get_col('dec')[0] 
        z_r = ran_selector.get_col('z')[0] 
        w_r = ran_selector.get_col('weight')[0] 
        print 'Original weights randoms:', w_r
        w_r = np.ones(len(ra_r))
        print 'Weights randoms, setting them to one:', w_r

        # Number of galaxies in each lens bin
        n_lens = np.array([len(ra_l[(z_l<zbins[i+1])&(z_l>zbins[i])]) for i in range(len(zbins)-1)])

        # Number of randoms in each redshift bin
        n_rand = np.array([len(ra_r[(z_r<zbins[i+1])&(z_r>zbins[i])]) for i in range(len(zbins)-1)])

        # Downsample randoms by this quantity: we want 10 times as many randoms as lenses per z-bin
        d = (n_lens*10.)/n_rand
        print 'd', d
        print 'n_lens', len(ra_l), n_lens
        print 'n_rand', len(ra_r), n_rand

        # Downsample
        np.random.seed(0)
        r = [np.random.rand(n_rand[i]) for i in range(len(n_rand))]
        ra_r = np.concatenate([ra_r[(z_r<zbins[i+1])&(z_r>zbins[i])][r[i] < d[i]] for i in range(len(zbins)-1)]) 
        dec_r = np.concatenate([dec_r[(z_r<zbins[i+1])&(z_r>zbins[i])][r[i] < d[i]] for i in range(len(zbins)-1)]) 
        z_r = np.concatenate([z_r[(z_r<zbins[i+1])&(z_r>zbins[i])][r[i] < d[i]] for i in range(len(zbins)-1)]) 

        print 'nrand after', len(ra_r)

        # Path to save the lens catalogs already jackknifed
        os.system('mkdir -p %s'%path_save)

        # Divides into jackknife regions
        jk_l = jk.jk(ra_l,dec_l,path_save)
        if type(jk_l) is int: 
                jk.jk_centers(ra_r,dec_r,path_save)
                jk_l = jk.jk(ra_l,dec_l,path_save)
        jk_r = jk.jk(ra_r,dec_r,path_save)

        print 'Number of lenses:', len(ra_l)

        c1 = fits.Column(name='RA', format='E', array=ra_l)
        c2 = fits.Column(name='DEC', format='E', array=dec_l)
        c3 = fits.Column(name='Z', format='E', array=z_l)
        c4 = fits.Column(name='ZERR', format='E', array=zerr_l)
        c5 = fits.Column(name='JK', format='K', array=jk_l)
        c6 = fits.Column(name='W', format='E', array=w_l)
        c7 = fits.Column(name='ID', format='K', array=ids)

        CC = [c1,c2,c3,c4,c5,c6,c7]
        #hdu = pf.new_table(CC, nrows=len(ra_l))
        t = fits.BinTableHDU.from_columns(CC)
        t.writeto('%s/lens.fits'%path_save, overwrite=True)

        c1 = fits.Column(name='RA', format='E', array=ra_r)
        c2 = fits.Column(name='DEC', format='E', array=dec_r)
        c3 = fits.Column(name='Z', format='E', array=z_r)
        c4 = fits.Column(name='JK', format='K', array=jk_r)
        CC = [c1,c2,c3,c4]
        t = fits.BinTableHDU.from_columns(CC)
        t.writeto('%s/random.fits'%path_save, overwrite=True)


print zbins
# Read yaml file that defines all the catalog selections used
params = yaml.load(open(paths['yaml']))
params['param_file'] = paths['yaml']

# Redmagic catalog
red_lens_selector, red_lens_calibrator = destest_functions.load_catalog(
    params, 'rm', None, params['lens_group'], params['lens_table'], params['lens_path'], return_calibrator=destest.NoCalib)

# Redmagic random catalog
red_ran_selector = destest_functions.load_catalog(
    params, 'ran', None, params['ran_group'], params['ran_table'], params['ran_path'])

# Mag lim catalog
alt_lens_selector, alt_lens_calibrator = destest_functions.load_catalog(
    params, 'maglim', None, params['alt_lens_group'], params['alt_lens_table'], params['alt_lens_path'], return_calibrator=destest.NoCalib)

# Mag lim random catalog
alt_ran_selector = destest_functions.load_catalog(
    params, 'ran', None, params['alt_ran_group'], params['alt_ran_table'], params['alt_ran_path'])

# DNF catalog for the alternative lens sample
pz_selector = destest_functions.load_catalog(
        params, 'pzdnf', None, params['dnf_group'], params['dnf_table'], params['dnf_path'], inherit=alt_lens_selector)

# Gold catalog for the alternative lens sample
gold_selector = destest_functions.load_catalog(
        params, 'gold', 'mcal', params['gold_group'], params['gold_table'], params['gold_path'], inherit=alt_lens_selector)


# First create the redmagic lens.fits and random.fits
if 'redmagic' in config['lens_v']:
        create_lens(red_lens_selector, red_ran_selector, None, None, zbins, paths['lenscats'], 'redmagic')

# Create lens.fits and random.fits for maglimit sample
if 'maglim' in config['lens_v']:
        create_lens(alt_lens_selector, alt_ran_selector, pz_selector, gold_selector, zbins, paths['lenscats'], 'maglim')
