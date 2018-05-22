import numpy as np                                                                               
import jk
from esutil.coords import randsphere
import sys
sys.path.append('../../destest/')                                                                
sys.path.append('../../kmeans_radec/')                                                           
import pyfits as pf
from info import zbins
import kmeans_radec
import os
import jk
import process_mice

path_mice = '../../../y1_shear_tests/sims/mice/'                                                  
files = [f for f in os.listdir(path_mice) if 'redmagic' in f]                                     
ra = np.zeros(0)
dec = np.zeros(0)
z = np.zeros(0)
for f in files:
	data = pf.getdata(path_mice + f)	
	ra = np.append(ra,data['ra'])
	dec = np.append(dec,data['dec'])
	z = np.append(z,data['ZREDMAGIC'])
    
np.random.seed(0)
r = np.random.rand(len(ra))
maskr = r<0.1
ra = ra[maskr]
dec = dec[maskr]
z = z[maskr]

# Load N(z)s from y1 redmagic data
path_cats = '/Volumes/Data/y1_shear_tests/cats/'
nzs = np.loadtxt(path_cats + 'redmagic/lens')
nzs = nzs.T
nz_zmean = nzs[0]
nz_zlims = np.zeros(len(nz_zmean) + 1)
nz_zlims[:-1] = nz_zmean - (nz_zmean[2] - nz_zmean[1])/2. 
nz_zlims[-1] = nz_zmean[-1] + (nz_zmean[2] - nz_zmean[1])/2.
nz_all = []
for l in range(0,5):
    nz_all.append(nzs[l+1])

nz_all = np.array(nz_all)    
zbins, nz_weights = process_mice.get_zbins_match_nz(zs = z_true, nz_zlims=nz_zlims, nzs=nz_all)





# Generating the randoms
ra_min = 0. 
dec_min = 0.
ra_max = 90.
dec_max = 90.
n_rand = len(ra)*2
ra_r, dec_r = randsphere(n_rand, ra_range=[ra_min,ra_max], dec_range=[dec_min, dec_max])

jk_l = jk.jk(ra,dec,path_mice)
if type(jk_l) is int: 
	jk.jk_centers(ra_r,dec_r,path_mice)
	jk_l = jk.jk(ra,dec,path_mice)
jk_r = jk.jk(ra_r,dec_r,path_mice)

c1 = pf.Column(name='RA', format='E', array=ra)
c2 = pf.Column(name='DEC', format='E', array=dec)
c3 = pf.Column(name='Z', format='E', array=z)
c4 = pf.Column(name='JK', format='K', array=jk_l)
CC = [c1,c2,c3,c4]
hdu = pf.new_table(CC, nrows=len(ra))
hdu.writeto('%s/lens.fits'%path_mice, clobber=True)

c1 = pf.Column(name='RA', format='E', array=ra_r)
c2 = pf.Column(name='DEC', format='E', array=dec_r)
c3 = pf.Column(name='JK', format='K', array=jk_r)
CC = [c1,c2,c3]
hdu = pf.new_table(CC, nrows=len(ra_r))
hdu.writeto('%s/random.fits'%path_mice, clobber=True)


