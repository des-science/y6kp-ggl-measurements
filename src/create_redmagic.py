import numpy as np
import sys
import yaml
sys.path.append('../../destest/')
sys.path.append('../../kmeans_radec/')
import destest
import pyfits as pf
from info import zbins, config, filename_mastercat, paths, filename_mastercat
import kmeans_radec
import os
import jk


mcal_file = paths['yaml'] + 'destest_mcal.yaml'
params_mcal = yaml.load(open(mcal_file))
params_mcal['param_file'] = mcal_file
params_mcal['filename'] = filename_mastercat
source_mcal = destest.H5Source(params_mcal)
source_selector = destest.Selector(params_mcal,source_mcal)
source_calibrator = destest.MetaCalib(params_mcal,source_selector)

gold_file = paths['yaml'] + 'destest_gold.yaml'
params_gold = yaml.load(open(gold_file))
params_gold['param_file'] = gold_file
params_gold['filename'] = filename_mastercat
source_gold = destest.H5Source(params_gold)
gold_selector = destest.Selector(params_gold,source_gold,inherit=source_selector)


lens_file = paths['yaml'] + 'destest_redmagic.yaml'
params_lens = yaml.load(open(lens_file))
params_lens['param_file'] = lens_file
params_lens['filename'] = filename_mastercat
params_lens['table'] = [config['redmagic_v']]
params_lens['select_path'] = 'index/redmagic/%s/select'%config['redmagic_v']
source_lens = destest.H5Source(params_lens)
lens_selector = destest.Selector(params_lens,source_lens)
lens_calibrator = destest.NoCalib(params_lens,lens_selector)

w_l = lens_calibrator.calibrate('x',weight_only=True)
print 'w_l', w_l
print sl;fdja
gmask = lens_calibrator.selector.get_match()
# ra, dec from Gold have better precision
# They have same ordering as the redmagic columns (not currently, there is a bug)
# When it is fixed used the two lines below:
# ra_lgold = gold_selector.source.read('ra')[0][gmask][lens_calibrator.selector.get_mask()[0]]
# dec_lgold = gold_selector.source.read('dec')[0][gmask][lens_calibrator.selector.get_mask()[0]]
# For the moment we use ra and dec from redmagic:
ra_l = lens_selector.get_col('ra')[0]
dec_l = lens_selector.get_col('dec')[0]

'''
print 'ra_lgold', ra_lgold
print 'ra_l', ra_l
print 'dec_lgold', dec_lgold
print 'dec_l', dec_l
print len(ra_l)

print (ra_l-ra_lgold).mean()
print (ra_l-ra_lgold).std()
print (dec_l-dec_lgold).mean()
print (dec_l-dec_lgold).std()
'''

z_l = lens_selector.get_col('zredmagic')[0]
zerr_l = lens_selector.get_col('zredmagic_e')[0]
ids = lens_selector.get_col('coadd_object_id')[0]
print len(ids)
assert len(ra_l)==len(ids), 'Something is wrong.' 

zbins = zbins['lims']
maskzl = ((z_l>zbins[0])&(z_l<zbins[-1]))
ra_l = ra_l[maskzl]
dec_l = dec_l[maskzl]
z_l = z_l[maskzl]
zerr_l = zerr_l[maskzl]
ids = ids[maskzl]
if w_l==1:
    w_l = [w_l]*len(ra_l)

param_file = paths['yaml'] + 'destest_random.yaml'
params_rmr = yaml.load(open(param_file))
params_rmr['param_file'] = param_file
params_rmr['filename'] = filename_mastercat
params_rmr['table'] = [config['redmagic_v']]
params_rmr['select_path'] = 'index/redmagic/%s/random_select'%config['redmagic_v']
source_rmr = destest.H5Source(params_rmr)
selector_rmr = destest.Selector(params_rmr,source_rmr)

ra_r = selector_rmr.get_col('ra')[0] 
dec_r = selector_rmr.get_col('dec')[0] 
z_r = selector_rmr.get_col('z')[0] 
w_r = selector_rmr.get_col('weight')[0] 
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
path = '../lens_cats/redmagic/%s/%s/njk_%d/'%(filename_mastercat[53:-3], params_lens['table'][0],config['njk'])
os.system('mkdir -p %s'%path)

# Divides into jackknife regions
jk_l = jk.jk(ra_l,dec_l,path)
if type(jk_l) is int: 
	jk.jk_centers(ra_r,dec_r,path)
	jk_l = jk.jk(ra_l,dec_l,path)
jk_r = jk.jk(ra_r,dec_r,path)

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
hdu.writeto('%s/lens.fits'%path, clobber=True)

c1 = pf.Column(name='RA', format='E', array=ra_r)
c2 = pf.Column(name='DEC', format='E', array=dec_r)
c3 = pf.Column(name='Z', format='E', array=z_r)
c4 = pf.Column(name='JK', format='K', array=jk_r)
CC = [c1,c2,c3,c4]
hdu = pf.new_table(CC, nrows=len(ra_r))
hdu.writeto('%s/random.fits'%path, clobber=True)

  
