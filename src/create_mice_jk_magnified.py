import numpy as np
import sys
import yaml
sys.path.append('../../kmeans_radec/')
import astropy.io.fits as fits
from info import zbins, config, paths
import kmeans_radec
import os
import jk
import ipdb

print config

def create_lens(lens, rand, zbins, path_save, sample):

        """
        Function that saves lens.fits and random.fits files, after splitting into JK regions.
        zbins: dictionary from info with the limits of the sample. Here we only use the low and high limits 
               to avoid selecting galaxies we will not use later.
        path_save: where the files are saved.
        sample: string, can be 'redmagic' or 'maglim'.
        """

        ra_l, dec_l, z_l, w_l = lens['RA'], lens['DEC'], lens['ZREDMAGIC'], lens['WEIGHT']
        zerr_l = lens['ZREDMAGIC_E']
        ids = lens['COADD_OBJECTS_ID']

        print 'Sample:', sample
        print 'Saving lens sample here:', path_save

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

        print 'Weights lenses:', w_l
        print 'Range of lens weights:', np.min(w_l), np.max(w_l)
        print 'Mean of lens weights:', np.mean(w_l)
        print 'Len zero lens weights:', len(w_l[w_l==0])
        
        #Load ra,dec, from random catalog 
        ra_r, dec_r, z_r, w_r = rand['RA'], rand['DEC'], rand['Z'], rand['WEIGHT']
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


        print 'd', d
        print 'n_lens', len(ra_l), n_lens
        print 'n_rand', len(ra_r), n_rand
        

        # Path to save the lens catalogs already jackknifed
        os.system('mkdir -p %s'%path_save)

        # Divides into jackknife regions
        jk_l = jk.jk(ra_l,dec_l,path_save)
        if type(jk_l) is int: 
                jk.jk_centers(ra_r,dec_r,path_save)
                jk_l = jk.jk(ra_l,dec_l,path_save)
        jk_r = jk.jk(ra_r,dec_r,path_save)

        c1 = fits.Column(name='RA', format='E', array=ra_l)
        c2 = fits.Column(name='DEC', format='E', array=dec_l)
        c3 = fits.Column(name='Z', format='E', array=z_l)
        c4 = fits.Column(name='ZERR', format='E', array=zerr_l)
        c5 = fits.Column(name='JK', format='K', array=jk_l)
        c6 = fits.Column(name='W', format='E', array=w_l)
        c7 = fits.Column(name='ID', format='K', array=ids)

        CC = [c1,c2,c3,c4,c5,c6,c7]
        t = fits.BinTableHDU.from_columns(CC)
        t.writeto('%s/lens.fits'%path_save, overwrite=True)

        c1 = fits.Column(name='RA', format='E', array=ra_r)
        c2 = fits.Column(name='DEC', format='E', array=dec_r)
        c3 = fits.Column(name='Z', format='E', array=z_r)
        c4 = fits.Column(name='JK', format='K', array=jk_r)
        CC = [c1,c2,c3,c4]
        t = fits.BinTableHDU.from_columns(CC)
        t.writeto('%s/random.fits'%path_save, overwrite=True)




l = fits.open(paths['mice_redmagic'])
lens = l[1].data

r = fits.open(paths['mice_randoms'])
rand = r[1].data

create_lens(lens, rand, zbins, paths['mice_jk'], sample='redmagic')


