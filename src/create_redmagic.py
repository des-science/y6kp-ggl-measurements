import numpy as np
import sys
import yaml
sys.path.append('../../destest/')
sys.path.append('../../kmeans_radec/')
import destest
import pyfits as pf
from info import zbins
import kmeans_radec
import os
import jk

param_file = 'destest_redmagic.yaml'
params_rm = yaml.load(open(param_file))
params_rm['param_file'] = param_file
source_rm = destest.H5Source(params_rm)
selector_rm = destest.Selector(params_rm,source_rm)

zbins = zbins['lims']
z_l = selector_rm.get_col('zredmagic')[0] 
print len(z_l)
maskzl = ((z_l>zbins[0])&(z_l<zbins[-1]))
z_l = z_l[maskzl]
ra_l = selector_rm.get_col('ra')[0][maskzl]
dec_l = selector_rm.get_col('dec')[0][maskzl] 
zerr_l = selector_rm.get_col('zredmagic_e')[0][maskzl] 
ids = selector_rm.get_col('coadd_object_id')[0][maskzl]
w_l = np.ones(len(ra_l))
print len(z_l)

param_file = 'destest_randoms.yaml'
params_rmr = yaml.load(open(param_file))
params_rmr['param_file'] = param_file
source_rmr = destest.H5Source(params_rmr)
selector_rmr = destest.Selector(params_rmr,source_rmr)

ra_r = selector_rmr.get_col('ra')[0] 
dec_r = selector_rmr.get_col('dec')[0] 
z_r = selector_rmr.get_col('z')[0] 
w_r = selector_rmr.get_col('weight')[0] 

# Number of galaxies in each lens bin
n_lens = np.array([len(ra_l[(z_l<zbins[i+1])&(z_l>zbins[i])]) for i in range(len(zbins)-1)])

# Number of randoms in each redshift bin
n_rand = np.array([len(ra_r[(z_r<zbins[i+1])&(z_r>zbins[i])]) for i in range(len(zbins)-1)])

# Downsample by this quantity 
d = 1/((n_rand/n_lens/10.))

# Downsample
np.random.seed(0)
r = [np.random.rand(n_rand[i]) for i in range(len(n_rand))]
ra_r = np.concatenate([ra_r[(z_r<zbins[i+1])&(z_r>zbins[i])][r[i] < d[i]] for i in range(len(zbins)-1)]) 
dec_r = np.concatenate([dec_r[(z_r<zbins[i+1])&(z_r>zbins[i])][r[i] < d[i]] for i in range(len(zbins)-1)]) 
z_r = np.concatenate([z_r[(z_r<zbins[i+1])&(z_r>zbins[i])][r[i] < d[i]] for i in range(len(zbins)-1)]) 

path = '../lens_cats/redmagic/%s'%params_rm['table'][0]
os.system('mkdir %s'%path)
with open("%s/version_name"%path, "w") as text_file:
    text_file.write("%s"%params_rm['filename'])

jk_l = jk.jk(ra_l,dec_l,path)
if type(jk_l) is int: 
	jk.jk_centers(ra_r,dec_r,path)
	jk_l = jk.jk(ra_l,dec_l,path)
jk_r = jk.jk(ra_r,dec_r,path)


print len(ra_l)
print sjldf
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

