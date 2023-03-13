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

        ra_l, dec_l, z_l = lens['ra_gal'], lens['dec_gal'], lens['z_cgal'] #z_cgal is ztrue, to bin in the same way as Ismael, so Georgios can compare the measurements
        ids = lens['coadd_objects_id']

        print 'Sample:', sample
        print 'Saving lens sample here:', path_save

        # Select galaxies that belong to some lens z-bin
        zbins = zbins['lims']
        print 'zbins:', zbins
        maskzl = ((z_l>zbins[0])&(z_l<zbins[-1]))
        ra_l = ra_l[maskzl]
        dec_l = dec_l[maskzl]
        z_l = z_l[maskzl]
        ids = ids[maskzl]
        
        #Load ra,dec, from random catalog 
        ra_r, dec_r, z_r = rand['ra'], rand['dec'], rand['z']

        # Number of galaxies in each lens bin
        n_lens = np.array([len(ra_l[(z_l<zbins[i+1])&(z_l>zbins[i])]) for i in range(len(zbins)-1)])

        # Number of randoms in each redshift bin
        n_rand = np.array([len(ra_r[(z_r<zbins[i+1])&(z_r>zbins[i])]) for i in range(len(zbins)-1)])

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
        c4 = fits.Column(name='JK', format='K', array=jk_l)
        c5 = fits.Column(name='ID', format='K', array=ids)

        CC = [c1,c2,c3,c4,c5]
        t = fits.BinTableHDU.from_columns(CC)
        t.writeto('%s/lens.fits'%path_save, overwrite=True)

        c1 = fits.Column(name='RA', format='E', array=ra_r)
        c2 = fits.Column(name='DEC', format='E', array=dec_r)
        c3 = fits.Column(name='Z', format='E', array=z_r)
        c4 = fits.Column(name='JK', format='K', array=jk_r)
        CC = [c1,c2,c3,c4]
        t = fits.BinTableHDU.from_columns(CC)
        t.writeto('%s/random.fits'%path_save, overwrite=True)

path_original_files = '../cats/mice/v_2/unmagnified/redmagic/'
lens_highdens = 'high_density_mice2_unmagnified.fits'
ran_highdens = 'mice2_desy3_v5_1_run_redmagic_highdens_randoms.fit'

# for the HOD project we only need the first three lens bins, so no need to use the highlum cats
#lens_highlum = 'high_lum_mice2_unmagnified.fits'
#ran_highlum = 'mice2_desy3_v5_1_run_redmagic_highlum_randoms.fit'

l = fits.open(path_original_files + lens_highdens)
lens = l[1].data

r = fits.open(path_original_files + ran_highdens)
rand = r[1].data

create_lens(lens, rand, zbins, paths['mice_jk'], sample='redmagic')


