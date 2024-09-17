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
        ra =  np.array(maglim['desy6kp//maglim/tomo_bin_{}'.format(zl_bin)]['RA'])
        dec = np.array(maglim['desy6kp//maglim/tomo_bin_{}'.format(zl_bin)]['DEC'])
        w   = np.array(maglim['desy6kp/maglim/tomo_bin_{}'.format(zl_bin)]['weight'])

        # weights = joblib.load('/global/cfs/cdirs/des/giannini/ggl/lss_weights_oct2022.pkl')
        # w = weights[zl_bin]

        maglim.close()
        del maglim
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
    
    
    def load_lens_Y6_maglim_old(self, path, zl_bin=None):
        """
        Loads lens galaxy data from file

        Options:
        - Y6 MagLim
        """
        # print ('HELLOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO')

        # read data for lenses
        hdul = fits.open(path+'/maglim++_bin{}.fits'.format(zl_bin))
        ra = hdul[1].data['RA']
        dec = hdul[1].data['DEC']
        # z = hdul[1].data['Z']
        w = hdul[1].data['W']

        hdul.close()
        
        return ra, dec, w
    

    def load_lens_Y3_maglim_orig(self, data_file, zl_lims=None):
        """
        Loads lens galaxy data from file

        Options:
        - Y3 MagLim
        """
        # print ('HELLOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO')

        # read data for lenses
        hdul = fits.open(data_file)
        ra = hdul[1].data['RA']
        dec = hdul[1].data['DEC']
        z = hdul[1].data['Z']
        w = hdul[1].data['W']

        # mask in redshift
        goodm = np.where( (z>zl_lims[0]) * (z<zl_lims[1]) )[0]
        
        hdul.close()
        
        return ra[goodm], dec[goodm], w[goodm]
    
    
    
    
    
    
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
    
    
    def load_randoms_Y6_maglim_old(self, path):
        """
        Loads random points data from file

        Options:
        - Y6 MagLim randoms
        """
        # read data for lenses
        hdul = fits.open(path+'/randoms++.fits')
        ra = hdul[1].data['RA']
        dec = hdul[1].data['DEC']

        hdul.close()
        
        return ra, dec
    

    def load_randoms_Y3_maglim_old(self, data_file, zl_lims=None):
        """
        Loads random points data from file

        Options:
        - Y3 MagLim randoms
        """
        # read data for lenses
        hdul = fits.open(data_file)
        ra = hdul[1].data['RA']
        dec = hdul[1].data['DEC']
        z = hdul[1].data['Z']

        # mask in redshift
        goodm = np.where( (z>zl_lims[0]) * (z<zl_lims[1]) )[0]
        
        hdul.close()
        
        return ra[goodm], dec[goodm]
    

    
    
    
    
    
    
    
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

    
    
    
    #Editted to load star data for psf residual tests
    
    def load_source_metadetect_old(self, path, zs_bin=None):
        """
        Loads Y6 source galaxy data
        
        """
            
        #resp = np.loadtxt(path+'/Response_bin{}.txt'.format(zs_bin))
        #file = pf.open(path+'/metadetect_bin{}.fits'.format(zs_bin))
        file = fits.open(path)
        
        ra_s = file[1].data['ra_s']
        dec_s = file[1].data['dec_s']
        e1_s = file[1].data['e1_s']
        e2_s = file[1].data['e2_s']
        #w_g = file[1].data['w_g']
        #R_g = resp[0]
        
#         mm['ra_s'] = ra_s 
#         mm['dec_s'] = dec_s
#         mm['e1_s'] = e1_s 
#         mm['e2_s'] = e2_s 
#         mm['g_cov_1_1'] = g_cov_1_1
#         mm['g_cov_2_2'] = g_cov_2_2
#         mm['w_g'] = w_g


#         # path = '/global/cfs/cdirs/des/giannini/ggl/metad_maglim_fidmask/'
#         path = par.out_main
#         if not os.path.exists(path):
#             os.makedirs(path)

#         meta_table = Table.from_pandas(mm)
#         meta_table.write(path+'/metadetect_bin{}.fits'.format(zs_bin), overwrite = True)

#         response = np.array([R_g, R11, R22])
#         np.savetxt(path+'/Response_bin{}.txt'.format(zs_bin), response)

        file.close()
    
        # return source['ra'] , source['dec'], source['e1'], source['e2'], R_g, w_g
        return ra_s, dec_s, e1_s, e2_s#, R_g, w_g

    
    
    
    
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
    
    
    
    
    
    
    
    

    def load_source_metadetect__(self, data_file, zs_bin=None):
        """
        Loads Y6 source galaxy data
        
        """
        print ('zs_bin', zs_bin)
        # read fits file metadetect tomographic bin
        met = fits.open(data_file)

        source = {}
        source['ra']  = np.array(met[1].data['ra_s'])
        source['dec'] = np.array(met[1].data['dec_s'])
        source['e1']  = np.array(met[1].data['e1_s'])
        source['e2']  = np.array(met[1].data['e2_s'])
        source['g_cov_1_1']  = np.array(met[1].data['g_cov_1_1'])
        source['g_cov_2_2']  = np.array(met[1].data['g_cov_2_2'])
        
        met.close()
        
        response = np.loadtxt('/global/cfs/cdirs/des/giannini/ggl/Response_bin{}.txt'.format(zs_bin))
        R_g = response[0]
        
        def get_shear_weights(cols, weights_file, weight_type):

            def _assign_loggrid(x, y, xmin, xmax, xsteps, ymin, ymax, ysteps):
                from math import log10
                # return x and y indices of data (x,y) on a log-spaced grid that runs from [xy]min to [xy]max in [xy]steps

                logstepx = log10(xmax/xmin)/xsteps
                logstepy = log10(ymax/ymin)/ysteps

                indexx = (np.log10(x/xmin)/logstepx).astype(int)
                indexy = (np.log10(y/ymin)/logstepy).astype(int)

                indexx = np.maximum(indexx,0)
                indexx = np.minimum(indexx, xsteps-1)
                indexy = np.maximum(indexy,0)
                indexy = np.minimum(indexy, ysteps-1)

                return indexx,indexy

            def _find_shear_weight(d, wgt_dict, snmin, snmax, sizemin, sizemax, steps):

                if wgt_dict is None:
                    weights = np.ones(len(d))
                    return weights

                shear_wgt = wgt_dict['weight']
                indexx, indexy = _assign_loggrid(d['s2n'], d['T_ratio'], snmin, snmax, steps, sizemin, sizemax, steps)
                weights = np.array([shear_wgt[x, y] for x, y in zip(indexx, indexy)])

                return weights

            #import pdb ; pdb.set_trace()
            if weight_type == 's2n_sizer':
                # pickle file that defines w(S/N, size)
                with open(weights_file, 'rb') as handle:
                    wgt_dict = pickle.load(handle)
                ## TO-DO: make snmin, snmax, sizemin, sizemax available in config file. 
                shear_wgt = _find_shear_weight(cols, wgt_dict, 10, 1000, 0.5, 5.0, 20)
            elif weight_type == 'shape_err':
                shear_wgt = 1/(0.17**2 + 0.5*(cols['g_cov_1_1'] + cols['g_cov_2_2']))

            shear_wgt[np.isnan(shear_wgt)] = 0.

            return shear_wgt
    
        weights_file = '/global/cfs/cdirs/des/myamamot/y6_shear_catalogs/Y6A2_METADETECT_V4/inverse_variance_weight_v3_s2n_10-300_Tratio_0.5-5.pickle'
        weight_type = 'shape_err'
        
        
        source['w'] = get_shear_weights(source, weights_file, weight_type)
        w_g = source['w']
        print ('wg', w_g[:10])

        return source['ra'] , source['dec'], source['e1'], source['e2'], R_g, w_g
    
    
    
    
 
  
    
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
                                 g1=(e1), 
                                 g2=(e2), 
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
            nn_lp = treecorr.NNCorrelation(nbins=self.par.ang_nbins, min_sep=theta_min, max_sep=theta_max, sep_units='arcmin', bin_slop=self.par.bin_slop, var_method='jackknife')
            nn_rp = treecorr.NNCorrelation(nbins=self.par.ang_nbins, min_sep=theta_min, max_sep=theta_max, sep_units='arcmin', bin_slop=self.par.bin_slop, var_method='jackknife')
            # catalog containing single point
            cat_p = treecorr.Catalog(ra=np.array([np.mean(ra_l)]), dec=np.array([np.mean(dec_l)]), ra_units=units, dec_units=units)
            # process catalogs and get objects 
            # containing the weights split into jackknife patches
            nn_lp.process(cat_l, cat_p, low_mem=low_mem)
            nn_rp.process(cat_r, cat_p, low_mem=low_mem)

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
                func = lambda corrs: ( (corrs[0].weight/np.sum(corrs[2].weight)) / (corrs[1].weight/np.sum(corrs[3].weight)) 
                                        * corrs[0].xi/corrs[0].Rg )
                corrs = [ng,rg,nn_lp,nn_rp]
            else:
                func = lambda corrs: corrs[0].xi / corrs[0].Rg
                corrs = [ng]
        cov_jk_gt = treecorr.estimate_multi_cov(corrs, 'jackknife', func=func)

        # get gammax covariance
        if use_randoms:
            func = lambda corrs: corrs[0].xi_im/corrs[0].Rg - corrs[1].xi_im/corrs[1].Rg
            corrs = [ng,rg]
        else:
            func = lambda corrs: corrs[0].xi_im/corrs[0].Rg
            corrs = [ng]
        cov_jk_gx = treecorr.estimate_multi_cov(corrs, 'jackknife', func=func)

        # get boost factor covariance
        if use_boosts:
            func = lambda corrs: (corrs[0].weight/np.sum(corrs[2].weight)) / (corrs[1].weight/np.sum(corrs[3].weight))
            corrs = [ng,rg,nn_lp,nn_rp]
            cov_jk_boost = treecorr.estimate_multi_cov(corrs, 'jackknife', func=func)
        else:
            cov_jk_boost = np.zeros((len(theta),len(theta)))

        # get covariance of randoms points
        corrs = [rg]
        func = lambda corrs: corrs[0].xi/corrs[0].Rg
        cov_jk_rand = treecorr.estimate_multi_cov(corrs, 'jackknife', func=func)

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