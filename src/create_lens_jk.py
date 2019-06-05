import numpy as np
import sys
import yaml
sys.path.append('../../destest/')
sys.path.append('../../kmeans_radec/')
import destest
import pyfits as pf
from info import zbins, config, paths
import destest_functions
import kmeans_radec
import os
import jk

# Read yaml file that defines all the catalog selections used
params = yaml.load(open(paths['yaml']))
params['param_file'] = paths['yaml']

# Redmagic catalog
lens_selector, lens_calibrator = destest_functions.load_catalog(
    params, None, params['lens_group'], params['lens_table'], params['lens_path'], return_calibrator=destest.NoCalib)

# Redmagic random catalog
ran_selector = destest_functions.load_catalog(
    params, None, params['ran_group'], params['ran_table'], params['ran_path'])

# Mag lim catalog
alt_lens_selector, alt_lens_calibrator = destest_functions.load_catalog(
    params, None, params['alt_lens_group'], params['alt_lens_table'], params['alt_lens_path'], return_calibrator=destest.NoCalib)

# Load ra,dec from redmagic catalog
ra_l  = lens_selector.get_col('ra')[0]
dec_l = lens_selector.get_col('dec')[0]
z_l = lens_selector.get_col('zredmagic')[0]
zerr_l = lens_selector.get_col('zredmagic_e')[0]
ids = lens_selector.get_col('coadd_object_id')[0]
w_l = lens_selector.get_col('weight')[0]
print w_l
assert len(ra_l)==len(ids), 'Something is wrong.' 

# Select galaxies that belong to some lens z-bin
zbins = zbins['lims']
maskzl = ((z_l>zbins[0])&(z_l<zbins[-1]))
ra_l = ra_l[maskzl]
dec_l = dec_l[maskzl]
z_l = z_l[maskzl]
zerr_l = zerr_l[maskzl]
ids = ids[maskzl]
w_l = w_l[maskzl]

#Load ra,dec, from random catalog of redmagic
ra_r = ran_selector.get_col('ra')[0] 
dec_r = ran_selector.get_col('dec')[0] 
z_r = ran_selector.get_col('z')[0] 
w_r = ran_selector.get_col('weight')[0] 
print 'w_r', w_r

# Number of galaxies in each lens bin
n_lens = np.array([len(ra_l[(z_l<zbins[i+1])&(z_l>zbins[i])]) for i in range(len(zbins)-1)])

# Number of randoms in each redshift bin
n_rand = np.array([len(ra_r[(z_r<zbins[i+1])&(z_r>zbins[i])]) for i in range(len(zbins)-1)])

# Downsample randoms by this quantity: we want 10 times as many randoms as lenses per z-bin
d = 1/((n_rand/n_lens/10.))

# Downsample
np.random.seed(0)
r = [np.random.rand(n_rand[i]) for i in range(len(n_rand))]
ra_r = np.concatenate([ra_r[(z_r<zbins[i+1])&(z_r>zbins[i])][r[i] < d[i]] for i in range(len(zbins)-1)]) 
dec_r = np.concatenate([dec_r[(z_r<zbins[i+1])&(z_r>zbins[i])][r[i] < d[i]] for i in range(len(zbins)-1)]) 
z_r = np.concatenate([z_r[(z_r<zbins[i+1])&(z_r>zbins[i])][r[i] < d[i]] for i in range(len(zbins)-1)]) 

# Path to save the lens catalogs already jackknifed
path_redmagic = paths['redmagic', config['redmagic_v']]
os.system('mkdir -p %s'%path_redmagic)

# Divides into jackknife regions
jk_l = jk.jk(ra_l,dec_l,path_redmagic)
if type(jk_l) is int: 
	jk.jk_centers(ra_r,dec_r,path_redmagic)
	jk_l = jk.jk(ra_l,dec_l,path_redmagic)
jk_r = jk.jk(ra_r,dec_r,path_redmagic)

print 'Number of lenses:', len(ra_l)
c1 = pf.Column(name='RA', format='E', array=ra_l)
c2 = pf.Column(name='DEC', format='E', array=dec_l)
c3 = pf.Column(name='Z', format='E', array=z_l)
c4 = pf.Column(name='ZERR', format='E', array=zerr_l)
c5 = pf.Column(name='JK', format='K', array=jk_l)
c6 = pf.Column(name='W', format='E', array=w_l)
c7 = pf.Column(name='ID', format='K', array=ids)

CC = [c1,c2,c3,c4,c5,c6,c7]
hdu = pf.new_table(CC, nrows=len(ra_l))
hdu.writeto('%s/lens.fits'%path_redmagic, clobber=True)

c1 = pf.Column(name='RA', format='E', array=ra_r)
c2 = pf.Column(name='DEC', format='E', array=dec_r)
c3 = pf.Column(name='Z', format='E', array=z_r)
c4 = pf.Column(name='JK', format='K', array=jk_r)
CC = [c1,c2,c3,c4]
hdu = pf.new_table(CC, nrows=len(ra_r))
hdu.writeto('%s/random.fits'%path_redmagic, clobber=True)

  
