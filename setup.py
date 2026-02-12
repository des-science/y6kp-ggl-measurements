"""
----------------
Helper functions
----------------
"""
import numpy as np
import os
import sys
import yaml
import treecorr
import h5py as h5
from astropy.io import fits
import healsparse
import joblib
import healpy as hp
from astropy.table import Table
import pandas as pd
import gc



class GGL_setup(object):

    def __init__(self, input_dir=None):

        " setup the parameters file "
        if input_dir is not None:
            if os.path.exists(input_dir):
                sys.path.append(input_dir)
            else:
                errmsg = '!!!Error: The input directory %s does not exist'%input_dir
                raise Exception(errmsg)
        else:
            errmsg = '!!!Error: Please provide the path to the input directory'
            raise Exception(errmsg)
        
        import params as par
        self.par = par
                    
        return

    
    
    
    def load_lens_Y6_maglim(self, path, zl_bin=None):
        """
        Loads lens galaxy data from file

        Options:
        - Y6 MagLim
        """
                         
        # read data for lenses
        maglim = h5.File(path, 'r')
        ra =  np.array(maglim['desy6kp/maglim/tomo_bin_{}'.format(zl_bin)]['RA'])
        dec = np.array(maglim['desy6kp/maglim/tomo_bin_{}'.format(zl_bin)]['DEC'])
        w   = np.array(maglim['desy6kp/maglim/tomo_bin_{}'.format(zl_bin)]['weight'])

        # weights = joblib.load('/global/cfs/cdirs/des/giannini/ggl/lss_weights_oct2022.pkl')
        # w = weights[zl_bin]

        maglim.close()
        del maglim
        gc.collect()
        
        return ra, dec, w
    
    
    def load_randoms_as_maglim(self, path, zl_bin=None):
        """
        Loads random points data from file

        Options:
        - Y6 lens randoms
        """
        # read data for lenses
                                                
        randoms = h5.File(path, 'r')
        ra =  np.array(randoms['desy6kp/ran_maglim/tomo_bin_0']['ra'])
        dec = np.array(randoms['desy6kp/ran_maglim/tomo_bin_0']['dec'])
        w =   np.ones(len(dec))
        
        randoms.close()
        del randoms
        gc.collect()
        
        return ra, dec, w
    

    
    def load_lens_Y6_maglim_all_bins(self, path, zl_bin=None):
        """
        Loads lens galaxy data from file

        Options:
        - Y6 MagLim
        """
                         
        # read data for lenses
        maglim = h5.File(path, 'r')
        ra = np.array([])
        dec = np.array([])
        for b in range(6):
            ra = np.append(ra, np.array(maglim['desy6kp//maglim/tomo_bin_{}'.format(b)]['RA']))
            print(len(np.array(maglim['desy6kp//maglim/tomo_bin_{}'.format(b)]['RA'])))
            dec = np.append(dec, np.array(maglim['desy6kp//maglim/tomo_bin_{}'.format(b)]['DEC']))
        weights = joblib.load('/global/cfs/cdirs/des/giannini/ggl/lss_weights_oct2022.pkl')
        w = weights[zl_bin]
        # w = np.ones(len(ra))
        # w = maglim_bin['W']
        
        print("ra length = " + str(len(ra)))
        maglim.close()
        del maglim
        gc.collect()
        
        return ra, dec, w
    

    
    
    def load_randoms_Y6(self, path, zl_bin=None):
        """
        Loads random points data from file

        Options:
        - Y6 lens randoms
        """
        # read data for lenses
                                                
        randoms = h5.File(path, 'r')
        ra =  np.array(randoms['desy6kp/ran_maglim/tomo_bin_{}'.format(zl_bin)]['ra'])
        dec = np.array(randoms['desy6kp/ran_maglim/tomo_bin_{}'.format(zl_bin)]['dec'])

        randoms.close()
        del randoms
        gc.collect()
        
        return ra, dec
    
    
    def load_alternative_randoms_Y6(self, path, zl_bin=None):
        """
        Loads random points data from file

        Options:
        - Y6 lens randoms
        """
        # read data for lenses
                                                
        randoms = h5.File(path, 'r')
        ra =  np.array(randoms['tomo_bin_{}'.format(zl_bin)]['ra'])
        dec = np.array(randoms['tomo_bin_{}'.format(zl_bin)]['dec'])

        randoms.close()
        del randoms
        gc.collect()
        
        return ra, dec
    
    

    
    
    
        #  BFD
    
    def load_source_bfd(self, path, zs_bin=None):
        """
        Loads Y6 source BFD galaxy data
        
        """
            
        bfd = h5.File(path, 'r')
        print ('Loads Y6 source BFD galaxy data')

        ra_s =  np.array(bfd['desy6kp/bfd']['tomo_bin_{}'.format(zs_bin)]['ra'])
        dec_s = np.array(bfd['desy6kp/bfd']['tomo_bin_{}'.format(zs_bin)]['dec'])
        e1_s =  np.array(bfd['desy6kp/bfd']['tomo_bin_{}'.format(zs_bin)]['e1'])
        e2_s =  np.array(bfd['desy6kp/bfd']['tomo_bin_{}'.format(zs_bin)]['e2'])
        w_g =   np.ones(len(ra_s))
        R_g = 1.   # has to stay like this, this is BFD!!
        
        bfd.close()
        del bfd
        gc.collect()
        
        return ra_s, dec_s, e1_s, e2_s, R_g, w_g
    
 
 
    def load_source_bfd_fits(self, path, zs_bin=None):
        """
        Loads Y6 source BFD galaxy data
        
        """
            
        bfd = pf.open(path)
        print ('Loads Y6 source BFD galaxy data')

        ra_s =  np.array(bfd[1].data['ra_s'])
        dec_s = np.array(bfd[1].data['dec_s'])
        e1_s =  np.array(bfd[1].data['e1_s'])
        e2_s =  np.array(bfd[1].data['e2_s'])
        w_g =   np.array(bfd[1].data['w'])
        # w_g =   np.ones(len(bfd[1].data['w']))
        R_g = 1.
        
        bfd.close()
        del bfd
        gc.collect()
        
        return ra_s, dec_s, e1_s, e2_s, R_g, w_g
    
        
        
        
        
        
        
        
        
        #   METADETECT
    def load_source_metadetect(self, path, r_path, zs_bin=None):
        """
        Loads Y6 source METADETECT galaxy data
        
        """
            
        mdet = h5.File(path, 'r')
        print ('Loads Y6 source METADETECT galaxy data')
        ra_s =  np.array(mdet['desy6kp/mdet/noshear']['tomo_bin_{}'.format(zs_bin)]['ra'])
        dec_s = np.array(mdet['desy6kp/mdet/noshear']['tomo_bin_{}'.format(zs_bin)]['dec'])
        e1_s =  np.array(mdet['desy6kp/mdet/noshear']['tomo_bin_{}'.format(zs_bin)]['gauss_g_1'])
        e2_s =  np.array(mdet['desy6kp/mdet/noshear']['tomo_bin_{}'.format(zs_bin)]['gauss_g_2'])
        w_g =   np.array(mdet['desy6kp/mdet/noshear']['tomo_bin_{}'.format(zs_bin)]['w'])
        # w_g =   np.ones(len(ra_s))

        resp = np.loadtxt(r_path)
        R_g = resp[0]

        mdet.close()
        del mdet
        gc.collect()
        
        return ra_s, dec_s, e1_s, e2_s, R_g, w_g


    
    def load_source_metadetect_unblinded(self, path, r_path, zs_bin=None):
        """
        Loads Y6 source METADETECT galaxy data
        
        """
        mdet = h5.File(path, 'r')
        print ('Loads Y6 source METADETECT galaxy data')
        ra_s =  np.array(mdet['noshear']['tomo_bin_{}'.format(zs_bin)]['ra'])
        dec_s = np.array(mdet['noshear']['tomo_bin_{}'.format(zs_bin)]['dec'])
        e1_s =  np.array(mdet['noshear']['tomo_bin_{}'.format(zs_bin)]['gauss_g_1'])
        e2_s =  np.array(mdet['noshear']['tomo_bin_{}'.format(zs_bin)]['gauss_g_2'])
        w_g =   np.array(mdet['noshear']['tomo_bin_{}'.format(zs_bin)]['w'])
        # w_g =   np.ones(len(ra_s))

        resp = np.loadtxt(r_path)
        R_g = resp[0]

        mdet.close()
        del mdet
        gc.collect()
        
        return ra_s, dec_s, e1_s, e2_s, R_g, w_g

    

    def load_PSF_shape(self, path, r_path, zs_bin=None):
        """
        Loads Y6 source galaxy data
        
        """
        file = fits.open('/global/cfs/cdirs/des/giannini/ggl/psf_res_shape.fits')
        #resp = np.loadtxt(path+'/Response_bin{}.txt'.format(zs_bin))
        #file = pf.open(path+'/metadetect_bin{}.fits'.format(zs_bin))
        # file = fits.open(path)
        
        ra_s = file[1].data['ra_s']
        dec_s = file[1].data['dec_s']
        e1_s = file[1].data['e1_s']
        e2_s = file[1].data['e2_s']
        w_g = np.ones(len(ra_s))
        R_g = 1.
        
        file.close()
        del file
        gc.collect()
    
        # return source['ra'] , source['dec'], source['e1'], source['e2'], R_g, w_g
        return ra_s, dec_s, e1_s, e2_s, R_g, w_g

    def load_PSF_size(self, path, r_path, zs_bin=None):
        """
        Loads Y6 source galaxy data
        
        """
        file = fits.open('/global/cfs/cdirs/des/giannini/ggl/psf_res_size.fits')
        #resp = np.loadtxt(path+'/Response_bin{}.txt'.format(zs_bin))
        #file = pf.open(path+'/metadetect_bin{}.fits'.format(zs_bin))
        # file = fits.open(path)
        
        ra_s = file[1].data['ra_s']
        dec_s = file[1].data['dec_s']
        e1_s = file[1].data['e1_s']
        e2_s = file[1].data['e2_s']
        w_g = np.ones(len(ra_s))
        R_g = 1.
        
        file.close()
        del file
        gc.collect()
    
        # return source['ra'] , source['dec'], source['e1'], source['e2'], R_g, w_g
        return ra_s, dec_s, e1_s, e2_s, R_g, w_g

    
    #Editted to load star data for psf residual tests    


    
    
    
    
    def load_mdet_with_sheared_ellipticities(self, path, zs_bin=None):
        """
        Loads Y6 METADETECT sheared ellipticities needed for scale-dependant response test
        
        """
            
        mdet = h5.File(path, 'r')
        print ('Loads Y6 source METADETECT sheared ellipticities')

        ra1p_s =  np.array(mdet['desy6kp/mdet/1p']['tomo_bin_{}'.format(zs_bin)]['ra'])
        dec1p_s = np.array(mdet['desy6kp/mdet/1p']['tomo_bin_{}'.format(zs_bin)]['dec'])
        w1p = np.array(mdet['desy6kp/mdet/1p']['tomo_bin_{}'.format(zs_bin)]['w'])
        
        ra1m_s =  np.array(mdet['desy6kp/mdet/1m']['tomo_bin_{}'.format(zs_bin)]['ra'])
        dec1m_s = np.array(mdet['desy6kp/mdet/1m']['tomo_bin_{}'.format(zs_bin)]['dec'])
        w1m = np.array(mdet['desy6kp/mdet/1m']['tomo_bin_{}'.format(zs_bin)]['w'])
        
        ra2p_s =  np.array(mdet['desy6kp/mdet/2p']['tomo_bin_{}'.format(zs_bin)]['ra'])
        dec2p_s = np.array(mdet['desy6kp/mdet/2p']['tomo_bin_{}'.format(zs_bin)]['dec'])
        w2p = np.array(mdet['desy6kp/mdet/2p']['tomo_bin_{}'.format(zs_bin)]['w'])
        
        ra2m_s =  np.array(mdet['desy6kp/mdet/2m']['tomo_bin_{}'.format(zs_bin)]['ra'])
        dec2m_s = np.array(mdet['desy6kp/mdet/2m']['tomo_bin_{}'.format(zs_bin)]['dec'])
        w2m = np.array(mdet['desy6kp/mdet/2m']['tomo_bin_{}'.format(zs_bin)]['w'])
        
        #for diagonal terms
        g1p_d = np.array(mdet['desy6kp/mdet/1p']['tomo_bin_{}'.format(zs_bin)]['gauss_g_1'])
        g1m_d = np.array(mdet['desy6kp/mdet/1m']['tomo_bin_{}'.format(zs_bin)]['gauss_g_1'])
        g2p_d = np.array(mdet['desy6kp/mdet/2p']['tomo_bin_{}'.format(zs_bin)]['gauss_g_2'])
        g2m_d = np.array(mdet['desy6kp/mdet/2m']['tomo_bin_{}'.format(zs_bin)]['gauss_g_2'])
        
        #for non-diagonal terms
        g1p_nd = np.array(mdet['desy6kp/mdet/1p']['tomo_bin_{}'.format(zs_bin)]['gauss_g_2'])
        g1m_nd = np.array(mdet['desy6kp/mdet/1m']['tomo_bin_{}'.format(zs_bin)]['gauss_g_2'])
        g2p_nd = np.array(mdet['desy6kp/mdet/2p']['tomo_bin_{}'.format(zs_bin)]['gauss_g_1'])
        g2m_nd = np.array(mdet['desy6kp/mdet/2m']['tomo_bin_{}'.format(zs_bin)]['gauss_g_1'])
        
        print("LENGTHS")
        print(len(g1p_d))
        print(len(g1m_d))
        print(len(g2p_d))
        print(len(g2m_d))
        print(len(w1p))

        mdet.close()
        del mdet
        gc.collect()
        
        return (ra1p_s, dec1p_s, ra1m_s, dec1m_s, ra2p_s, dec2p_s, ra2m_s, dec2m_s,
                w1p, w1m, w2p, w2m,
                g1p_d, g1m_d, g2p_d, g2m_d, g1p_nd, g1m_nd, g2p_nd, g2m_nd)
    
    
    
    
    
    
 
  
    
    def boost_factor_calculate(self, sum_w_l, sum_w_r, w_LS, w_RS):
        """
        Calculates boost factors

        input
        -----
        - w_l: weights associated with lenses
        - w_r: weights associated with randoms
        - w_LS: weights associated with lens-source pairs
        - w_RS: weights associated with random-source pairs

        output
        ------
        the boost factors as a function of theta
        """
        # boost factor calculation
        boosts = (sum_w_r/sum_w_l) * (w_LS/w_RS)

        return boosts


        
    
    def get_gammat(self, ra_l, dec_l, ra_s, dec_s, ra_rand=None, dec_rand=None, params=None, 
                   units='deg', sep_units='arcmin', low_mem=False, weights=None, 
                   use_randoms=False, use_boosts=False):
        """
        Calculates gamma_t as a function of theta using Treecorr

        input
        -----
        - ra_l, dec_l: ra and dec of lens galaxies
        - ra_rand, dec_rand: ra and dec of random points
        - ra_s, dec_s: ra and dec of source galaxies

        output
        ------
        theta, correlation functions, covariances, and extra information
        """
        # minimum and maximum angular separation
        # theta_min, theta_max = self.par.theta_lims

        # minimum and maximum angular separation
        nbins = 30 
        bin_edges = np.geomspace(2.5, 2500, nbins+1)
        bin_edges = bin_edges[:-4]
        nbins = 26
        bin_edges = np.geomspace(2.5, bin_edges[-1], nbins+1)
        print(bin_edges)
        theta_min = bin_edges[0]
        theta_max = bin_edges[-1]

        # ellipticity (e1,e2) and R_gamma
        e1, e2, Rg, wg = params

        # count-shear two-point correlation function(i.e. galaxy-galaxy lensing)
        print('Rg = ', Rg)
        print('Num lenses = ', len(ra_l))
        print('Num sources = ', len(ra_s))
        print('Average e1 = ', np.average(e1, weights=wg))
        print('Average e2 = ', np.average(e2, weights=wg))

        # generate lens catalogs to correlate and process them
        ng = treecorr.NGCorrelation(nbins=self.par.ang_nbins, min_sep=theta_min, max_sep=theta_max, sep_units=sep_units, bin_slop=self.par.bin_slop)
        cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l, ra_units=units, dec_units=units, w=weights)
        cat_s = treecorr.Catalog(ra=ra_s, dec=dec_s, ra_units=units, dec_units=units,
                                 g1=(e1), 
                                 g2=(e2), 
                                 w=wg)
        ng.process(cat_l, cat_s, low_mem=low_mem)

        # get the theta, gammat and covariance
        theta = np.exp(ng.logr)
        gamma_t = ng.xi / Rg

        # generate randoms catalogs to correlate and process them
        if use_randoms or use_boosts:
            rg = treecorr.NGCorrelation(nbins=self.par.ang_nbins, min_sep=theta_min, max_sep=theta_max, sep_units=sep_units, bin_slop=self.par.bin_slop)
            cat_r = treecorr.Catalog(ra=ra_rand, dec=dec_rand, ra_units=units, dec_units=units)
            rg.process(cat_r, cat_s, low_mem=low_mem)

        # get boost factors
        sum_w_l = np.sum(weights)
        if ra_rand is not None:
            sum_w_r = len(ra_rand)
        else:
            sum_w_r = 0
            
        if use_boosts:
            boost = self.boost_factor_calculate(sum_w_l, sum_w_r, ng.weight, rg.weight)
            gamma_t *= boost
        else:
            boost = np.ones(len(theta))

        if use_randoms:
            gammat_rand = rg.xi/Rg
            gamma_t -= gammat_rand
            xi_im_rand = rg.xi_im
            xi_npairs_rand = rg.npairs
            xi_weight_rand = rg.weight
        else:
            gammat_rand = np.zeros(len(gamma_t))
            xi_im_rand = np.zeros(len(gamma_t))
            xi_npairs_rand = np.zeros(len(gamma_t))
            xi_weight_rand = np.zeros(len(gamma_t))

        return (theta, gamma_t, ng.xi/Rg, gammat_rand, 
                ng.varxi, 
                ng.xi_im, xi_im_rand, ng.npairs, xi_npairs_rand, ng.weight, xi_weight_rand, 
                Rg, sum_w_l, sum_w_r, boost)


    
    def get_gammat_and_covariance(self, ra_l, dec_l, ra_s, dec_s, ra_rand=None, dec_rand=None, params=None, 
                                  units='deg', sep_units='arcmin', low_mem=False, weights=None, 
                                  use_randoms=False, use_boosts=False):
        """
        Calculates gamma_t as a function of theta and its Jackknife covariance using only Treecorr

        input
        -----
        - ra_l, dec_l: ra and dec of lens galaxies
        - ra_rand, dec_rand: ra and dec of random points
        - ra_s, dec_s: ra and dec of source galaxies

        output
        ------
        theta, correlation functions, covariances, and extra information
        """
        # minimum and maximum angular separation
        # theta_min, theta_max = self.par.theta_lims

        # minimum and maximum angular separation
        nbins = 30 
        bin_edges = np.geomspace(2.5, 2500, nbins+1)
        bin_edges = bin_edges[:-4]
        nbins = 26
        bin_edges = np.geomspace(2.5, bin_edges[-1], nbins+1)
        print(bin_edges)
        theta_min = bin_edges[0]
        theta_max = bin_edges[-1]
        # ellipticity (e1,e2), R_gamma and weights
        e1, e2, Rg, wg = params

        # count-shear two-point correlation function(i.e. galaxy-galaxy lensing)
        print('Rg = ', Rg)
        print('Num lenses = ', len(ra_l))
        print('Num sources = ', len(ra_s))
        print('Average e1 = ', np.average(e1, weights=wg))
        print('Average e2 = ', np.average(e2, weights=wg))

        # generate lens catalogs to correlate and process them
        ng = treecorr.NGCorrelation(nbins=self.par.ang_nbins, min_sep=theta_min, max_sep=theta_max, sep_units=sep_units, bin_slop=self.par.bin_slop, var_method='jackknife')
        print ('ng done')
        
        if os.path.isfile('jk_centers'):
            cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l, ra_units=units, dec_units=units, w=weights, patch_centers='jk_centers')
            print ('prepared catalog - loaded jk centers')
        else:
            cat_r = treecorr.Catalog(ra=ra_rand, dec=dec_rand, ra_units=units, dec_units=units, npatch=self.par.n_jck)
            print ('prepared random cat')
            cat_r.write_patch_centers('jk_centers')
            print ('wrote jk centers')
            cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l, ra_units=units, dec_units=units, w=weights, patch_centers='jk_centers')
            print ('done lens cat')


#         cat_r = treecorr.Catalog(ra=ra_rand, dec=dec_rand, ra_units=units, dec_units=units, npatch=self.par.n_jck)
#         print ('prepared random cat')
#         # make it save a txt file at the same time. 
#         np.savetxt('/global/homes/g/giannini/gglensing/y6kp-ggl-measurements/source/test.txt', np.zeros(1))
        
#         cat_r.write_patch_centers('jk_centers')
#         print ('wrote jk centers')
#         cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l, ra_units=units, dec_units=units, w=weights, patch_centers='jk_centers')
#         print ('done lens cat')


        # wg = np.ones(len(dec_s))
        cat_s = treecorr.Catalog(ra=ra_s, dec=dec_s, ra_units=units, dec_units=units,
                                 g1=(e1-np.average(e1, weights=wg)), 
                                 g2=(e2-np.average(e2, weights=wg)), 
                                 w=wg, patch_centers='jk_centers')
        print ('done source cat')
        ng.process(cat_l, cat_s, low_mem=low_mem)

        # get theta, gammat
        theta = np.exp(ng.logr)
        gamma_t = ng.xi/Rg              #For response test, compute Rg for each theta instead
        gammat_tot = np.copy(gamma_t)

        # get imaginary part of xi and gamma_x
        gamma_x = ng.xi_im/Rg
        gammax_tot = np.copy(gamma_x)

        # generate randoms catalogs to correlate and process them
        if use_randoms or use_boosts:
            rg = treecorr.NGCorrelation(nbins=self.par.ang_nbins, min_sep=theta_min, max_sep=theta_max, sep_units=sep_units, bin_slop=self.par.bin_slop, var_method='jackknife')
            
            cat_r = treecorr.Catalog(ra=ra_rand, dec=dec_rand, ra_units=units, dec_units=units, patch_centers='jk_centers')
            rg.process(cat_r, cat_s, low_mem=low_mem)

        # boost factors for gammat
        sum_w_l = np.sum(weights)
        if ra_rand is not None:
            sum_w_r = len(ra_rand)
        else:
            sum_w_r = 0
        if use_boosts:
            boost = self.boost_factor_calculate(sum_w_l, sum_w_r, ng.weight, rg.weight)
        else:
            boost = np.ones(len(theta))
        gammat_tot *= boost

        # generate fake treecorr correlation objects for lenses and randoms 
        # that hold the weights for the boost factor covariance calculations
        if use_boosts:
            # initialize NN correlations for single point cross lenses and randoms
            nn_lp = treecorr.NNCorrelation(nbins=1, min_sep=theta_min, max_sep=6000, sep_units='arcmin', bin_slop=self.par.bin_slop, var_method='jackknife')
            nn_rp = treecorr.NNCorrelation(nbins=1, min_sep=theta_min, max_sep=6000, sep_units='arcmin', bin_slop=self.par.bin_slop, var_method='jackknife')
            # catalog containing single point
            cat_p = treecorr.Catalog(ra=np.array([np.mean(ra_l)]), dec=np.array([np.mean(dec_l)]), ra_units=units, dec_units=units)
            # process catalogs and get objects 
            # containing the weights split into jackknife patches
            nn_lp.process(cat_l, cat_p, low_mem=low_mem)
            nn_rp.process(cat_r, cat_p, low_mem=low_mem)
            nn_lp.write('/global/cfs/cdirs/des/giannini/ggl/nn_lp_10Jun2025.fits', write_patch_results=True)
            nn_lp.write('/global/cfs/cdirs/des/giannini/ggl/nn_rp_10Jun2025.fits', write_patch_results=True)
            
        # update correlations with responses to use in Jackknife mode
        ng.Rg = Rg*np.ones(len(theta))

        # random-point subtraction for gamma_t and gamma_x
        if use_randoms:
            #correct for response
            gammat_rand = rg.xi/Rg
            gammax_rand = rg.xi_im/Rg
            xi_im_rand = rg.xi_im
            xi_npairs_rand = rg.npairs
            xi_weight_rand = rg.weight
            rg.Rg = Rg*np.ones(len(theta))
        else:
            gammat_rand = np.zeros(len(theta))
            gammax_rand = np.zeros(len(theta))
            xi_im_rand = np.zeros(len(theta))
            xi_npairs_rand = np.zeros(len(theta))
            xi_weight_rand = np.zeros(len(theta))
        gammat_tot -= gammat_rand
        gammax_tot -= gammax_rand

        # get gamma_t gammat covariance
        if use_randoms:
            if use_boosts:
                func = lambda corrs: ( (corrs[0].weight/np.sum(corrs[2].weight)) / (corrs[1].weight/np.sum(corrs[3].weight)) * corrs[0].xi/corrs[0].Rg - corrs[1].xi/corrs[1].Rg )
                corrs = [ng,rg,nn_lp,nn_rp]
            else:
                func = lambda corrs: ( corrs[0].xi/corrs[0].Rg - corrs[1].xi/corrs[1].Rg )
                corrs = [ng,rg]
        else:
            if use_boosts:
                func = lambda corrs: ( (corrs[0].weight/np.sum(corrs[2].weight)) / (corrs[1].weight/np.sum(corrs[3].weight))  * corrs[0].xi/corrs[0].Rg )
                corrs = [ng,rg,nn_lp,nn_rp]
            else:
                func = lambda corrs: corrs[0].xi / corrs[0].Rg
                corrs = [ng]
        cov_jk_gt = treecorr.estimate_multi_cov(corrs, 'jackknife', func=func, cross_patch_weight='match')

        # get gammax covariance
        if use_randoms:
            func = lambda corrs: corrs[0].xi_im/corrs[0].Rg - corrs[1].xi_im/corrs[1].Rg
            corrs = [ng,rg]
        else:
            func = lambda corrs: corrs[0].xi_im/corrs[0].Rg
            corrs = [ng]
        cov_jk_gx = treecorr.estimate_multi_cov(corrs, 'jackknife', func=func, cross_patch_weight='match')

        # get boost factor covariance
        if use_boosts:
            func = lambda corrs: (corrs[0].weight/np.sum(corrs[2].weight)) / (corrs[1].weight/np.sum(corrs[3].weight))
            corrs = [ng,rg,nn_lp,nn_rp]
            cov_jk_boost = treecorr.estimate_multi_cov(corrs, 'jackknife', func=func, cross_patch_weight='match')
        else:
            cov_jk_boost = np.zeros((len(theta),len(theta)))

        # get covariance of randoms points
        if use_randoms:
            corrs = [rg]
            func = lambda corrs: corrs[0].xi/corrs[0].Rg
            cov_jk_rand = treecorr.estimate_multi_cov(corrs, 'jackknife', func=func, cross_patch_weight='match')
        else:
            cov_jk_rand = np.zeros((len(theta),len(theta)))

        return (theta, gamma_t, gammat_tot, gammat_rand, gamma_x, gammax_tot, gammax_rand, 
                cov_jk_gt, ng.varxi, cov_jk_boost, cov_jk_gx, cov_jk_rand,
                ng.xi_im, xi_im_rand, ng.npairs, xi_npairs_rand, ng.weight, xi_weight_rand, 
                Rg, sum_w_l, sum_w_r, boost)
        
    
    def NK(self, ra_l, dec_l, ra_s, dec_s, ra_rand=None, dec_rand=None, params=None, 
                units='deg', sep_units='arcmin', low_mem=False, weights=None, 
                use_randoms=False, use_boosts=False):
        
        # minimum and maximum angular separation
        theta_min, theta_max = self.par.theta_lims

        # ellipticity for use in NK correlations
        K, wg = params
        
        print("K len")
        print(len(K))
        print("ra_s len")
        print(len(ra_s))

        # generate lens catalogs to correlate and process them
        nk = treecorr.NKCorrelation(nbins=self.par.ang_nbins, min_sep=theta_min, max_sep=theta_max, sep_units=sep_units, bin_slop=self.par.bin_slop, var_method='jackknife')
        print ('nk done')
        
        if os.path.isfile('jk_centers'):
            cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l, ra_units=units, dec_units=units, w=weights, patch_centers='jk_centers')
            print ('prepared catalog - loaded jk centers')
        else:
            cat_r = treecorr.Catalog(ra=ra_rand, dec=dec_rand, ra_units=units, dec_units=units, npatch=self.par.n_jck)
            print ('prepared random cat')
            cat_r.write_patch_centers('jk_centers')
            print ('wrote jk centers')
            cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l, ra_units=units, dec_units=units, w=weights, patch_centers='jk_centers')
            print ('done lens cat')
            
        cat_s = treecorr.Catalog(ra=ra_s, dec=dec_s, ra_units=units, dec_units=units,
                                 k=K, w=wg, patch_centers='jk_centers')
        print ('done source cat')
        nk.process(cat_l, cat_s, low_mem=low_mem)

        theta = np.exp(nk.logr)
        K_ang = nk.xi
        
        print("K_ang")
        print(K_ang)
        
        return theta, K_ang