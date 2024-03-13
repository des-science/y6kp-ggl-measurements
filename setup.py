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
import healsparse as hsp
import joblib
import healpy as hp
from astropy.table import Table
import pandas as pd
import gc
import pickle


class GGL_setup(object):

    def __init__(self, input_dir=None):

        " Load the parameters/settings file "
        if input_dir is not None:
            if os.path.exists(input_dir):
                sys.path.append(input_dir)
            else:
                errmsg = '!!!Error: The input directory {!s} does not exist'.format(input_dir)
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
                   
        with h5.File(path, 'r') as maglim:
            ra =  np.array(maglim['desy6kp/maglim/tomo_bin_{}'.format(zl_bin)]['RA'])
            dec = np.array(maglim['desy6kp/maglim/tomo_bin_{}'.format(zl_bin)]['DEC'])
            
        weights = joblib.load('/global/cfs/cdirs/des/giannini/ggl/lss_weights_oct2022.pkl')
        w = weights[zl_bin]
        # w = np.ones(len(ra))

        gc.collect()
        
        return ra, dec, w
    

    def load_lens_Y3_maglim(self, data_file, zl_lims=None):
        """
        Loads lens galaxy data from file

        Options:
        - Y3 MagLim
        """
        
        with fits.open(data_file) as hdul:
            ra = hdul[1].data['RA']
            dec = hdul[1].data['DEC']
            z = hdul[1].data['Z']
            w = hdul[1].data['W']

        # mask in redshift
        goodm = np.where( (z>zl_lims[0]) * (z<zl_lims[1]) )[0]
        
        return ra[goodm], dec[goodm], w[goodm]
    
    
    def load_randoms_Y6(self, path, zl_bin=None):
        """
        Loads random points data from file

        Options:
        - Y6 lens randoms
        """
                                                
        with h5.File(path, 'r') as randoms:
            ra =  np.array(randoms['desy6kp/ran/tomo_bin_{}'.format(zl_bin)]['RA'])
            dec = np.array(randoms['desy6kp/ran/tomo_bin_{}'.format(zl_bin)]['DEC'])

        gc.collect()
        
        return ra, dec
    

    def load_randoms_Y3_maglim(self, data_file, zl_lims=None):
        """
        Loads random points data from file

        Options:
        - Y3 MagLim randoms
        """
        
        with fits.open(data_file) as hdul:
            ra = hdul[1].data['RA']
            dec = hdul[1].data['DEC']
            z = hdul[1].data['Z']

        # mask in redshift
        goodm = np.where( (z>zl_lims[0]) * (z<zl_lims[1]) )[0]
        
        return ra[goodm], dec[goodm]
    
    
    def load_source_bfd(self, path, path_binning, path_mask, zs_bin=None):
        """
        Loads Y6 source BFD galaxy data
        
        """
        # Temporary: : the bins assignment in the mastercat is wrong and e1,e2 in the bfd cat are wrong
        #with h5.File(path, 'r') as bfd:
        #    ra =  np.array(bfd['desy6kp/bfd/tomo_bin_{}'.format(zs_bin)]['ra'])
        #    dec = np.array(bfd['desy6kp/bfd/tomo_bin_{}'.format(zs_bin)]['dec'])
        #    e1 =  np.array(bfd['desy6kp/bfd/tomo_bin_{}'.format(zs_bin)]['e1'])
        #    e2 =  np.array(bfd['desy6kp/bfd/tomo_bin_{}'.format(zs_bin)]['e2'])
        #    P =  np.array(bfd['desy6kp/bfd/tomo_bin_{}'.format(zs_bin)]['P'])
        #    Q0 =  np.array(bfd['desy6kp/bfd/tomo_bin_{}'.format(zs_bin)]['Q0'])
        #    Q1 =  np.array(bfd['desy6kp/bfd/tomo_bin_{}'.format(zs_bin)]['Q1'])
        #    R00 =  np.array(bfd['desy6kp/bfd/tomo_bin_{}'.format(zs_bin)]['R00'])
        #    R01 =  np.array(bfd['desy6kp/bfd/tomo_bin_{}'.format(zs_bin)]['R01'])
        #    R11 =  np.array(bfd['desy6kp/bfd/tomo_bin_{}'.format(zs_bin)]['R11'])
        
        with open(path_binning, 'rb') as file:
            tomo_bin = pickle.load(file)
            
        with fits.open(path) as bfd:
            ra =  np.array(bfd[1].data['ra'])[np.where(tomo_bin==zs_bin)[0]]
            dec = np.array(bfd[1].data['dec'])[np.where(tomo_bin==zs_bin)[0]]
            
        print(len(ra))
        mask = hsp.HealSparseMap.read(path_mask)
        mask_bool = mask.get_values_pos(ra, dec)
    
        with fits.open(path) as bfd:
            ra =  ra[mask_bool]
            dec = dec[mask_bool]
            P = np.array(bfd[1].data['P'])[np.where(tomo_bin==zs_bin)[0]][mask_bool]
            Q0 = np.array(bfd[1].data['Q0'])[np.where(tomo_bin==zs_bin)[0]][mask_bool]
            Q1 = np.array(bfd[1].data['Q1'])[np.where(tomo_bin==zs_bin)[0]][mask_bool]
            R00 = np.array(bfd[1].data['R00'])[np.where(tomo_bin==zs_bin)[0]][mask_bool]
            R01 = np.array(bfd[1].data['R01'])[np.where(tomo_bin==zs_bin)[0]][mask_bool]
            R11 = np.array(bfd[1].data['R11'])[np.where(tomo_bin==zs_bin)[0]][mask_bool]
            
        pqr = (np.vstack(P, Q0, Q1, R00, R01, R11).T).astype(np.float64)
        logpqr = logPQR(pqr)

        Q0 = logpqr[:,1]
        Q1 = logpqr[:,2]
        R00 = logpqr[:,3]
        R01 = logpqr[:,4]
        R11 = logpqr[:,5]
        
        e1, e2 = approx_es_bfd(logpqr)
        
        w_g =   np.ones(len(ra))
             
        gc.collect()
        
        return ra, dec, e1, e2, P, Q0, Q1, R00, R01, R11, w_g
    
    
    def logPQR(pqr):
        """
        Convert an Nx6 array representing quadratic Taylor 
        expansion of p w.r.t. shear into Taylor expansion of log(P).
        Any input rows with p<=0 are put as zeros in output
        """
        out = np.zeros_like(pqr)
        p = pqr[:,0]
        good = p>0.
        q = pqr[:,1:3] / p[:,np.newaxis]
        r = pqr[:,3:] / p[:,np.newaxis]
        out[:,0] = np.log(p)
        out[:,1:3] = q
        r[:,0] -= q[:,0]*q[:,0]
        r[:,1] -= q[:,1]*q[:,0]
        r[:,2] -= q[:,1]*q[:,1]
        out[:,3:] = r
        out = np.where(good[:,np.newaxis],out,0.)
        return out
    
    
    def approx_es_bfd(logpqr):
        Q = logpqr[:,1:3]
        Rs = np.sum(logpqr[:,3:], axis=0, dtype=float)
        R = np.array([[Rs[0], Rs[1]], [Rs[1], Rs[2]]])
        invR = np.linalg.inv(R)
        Q *= (len(Q[:,0]))
        es = np.matmul(invR, Q.T)  
        return es[0], es[1]
    
 
    def load_source_metadetect(self, path, r_path, zs_bin=None):
        """
        Loads Y6 source METADETECT galaxy data
        
        """
            
        with h5.File(path, 'r') as mdet:
            ra_s =  np.array(mdet['desy6kp/mdet/noshear']['tomo_bin_{}'.format(zs_bin)]['ra'])
            dec_s = np.array(mdet['desy6kp/mdet/noshear']['tomo_bin_{}'.format(zs_bin)]['dec'])
            e1_s =  np.array(mdet['desy6kp/mdet/noshear']['tomo_bin_{}'.format(zs_bin)]['gauss_g_1'])
            e2_s =  np.array(mdet['desy6kp/mdet/noshear']['tomo_bin_{}'.format(zs_bin)]['gauss_g_2'])
            w_g =   np.array(mdet['desy6kp/mdet/noshear']['tomo_bin_{}'.format(zs_bin)]['w'])
        
        resp = np.loadtxt(r_path)
        R_g = resp[0]

        gc.collect()
        
        return ra_s, dec_s, e1_s, e2_s, R_g, w_g
    

#     def load_source_metadetect__(self, data_file, zs_bin=None):
#         """
#         Loads Y6 source galaxy data
        
#         """
#         print ('zs_bin', zs_bin)
#         # read fits file metadetect tomographic bin
#         met = fits.open(data_file)

#         source = {}
#         source['ra']  = np.array(met[1].data['ra_s'])
#         source['dec'] = np.array(met[1].data['dec_s'])
#         source['e1']  = np.array(met[1].data['e1_s'])
#         source['e2']  = np.array(met[1].data['e2_s'])
#         source['g_cov_1_1']  = np.array(met[1].data['g_cov_1_1'])
#         source['g_cov_2_2']  = np.array(met[1].data['g_cov_2_2'])
        
#         met.close()
        
#         response = np.loadtxt('/global/cfs/cdirs/des/giannini/ggl/Response_bin{}.txt'.format(zs_bin))
#         R_g = response[0]
        
#         def get_shear_weights(cols, weights_file, weight_type):

#             def _assign_loggrid(x, y, xmin, xmax, xsteps, ymin, ymax, ysteps):
#                 from math import log10
#                 # return x and y indices of data (x,y) on a log-spaced grid that runs from [xy]min to [xy]max in [xy]steps

#                 logstepx = log10(xmax/xmin)/xsteps
#                 logstepy = log10(ymax/ymin)/ysteps

#                 indexx = (np.log10(x/xmin)/logstepx).astype(int)
#                 indexy = (np.log10(y/ymin)/logstepy).astype(int)

#                 indexx = np.maximum(indexx,0)
#                 indexx = np.minimum(indexx, xsteps-1)
#                 indexy = np.maximum(indexy,0)
#                 indexy = np.minimum(indexy, ysteps-1)

#                 return indexx,indexy

#             def _find_shear_weight(d, wgt_dict, snmin, snmax, sizemin, sizemax, steps):

#                 if wgt_dict is None:
#                     weights = np.ones(len(d))
#                     return weights

#                 shear_wgt = wgt_dict['weight']
#                 indexx, indexy = _assign_loggrid(d['s2n'], d['T_ratio'], snmin, snmax, steps, sizemin, sizemax, steps)
#                 weights = np.array([shear_wgt[x, y] for x, y in zip(indexx, indexy)])

#                 return weights

#             #import pdb ; pdb.set_trace()
#             if weight_type == 's2n_sizer':
#                 # pickle file that defines w(S/N, size)
#                 with open(weights_file, 'rb') as handle:
#                     wgt_dict = pickle.load(handle)
#                 ## TO-DO: make snmin, snmax, sizemin, sizemax available in config file. 
#                 shear_wgt = _find_shear_weight(cols, wgt_dict, 10, 1000, 0.5, 5.0, 20)
#             elif weight_type == 'shape_err':
#                 shear_wgt = 1/(0.17**2 + 0.5*(cols['g_cov_1_1'] + cols['g_cov_2_2']))

#             shear_wgt[np.isnan(shear_wgt)] = 0.

#             return shear_wgt
    
#         weights_file = '/global/cfs/cdirs/des/myamamot/y6_shear_catalogs/Y6A2_METADETECT_V4/inverse_variance_weight_v3_s2n_10-300_Tratio_0.5-5.pickle'
#         weight_type = 'shape_err'
        
        
#         source['w'] = get_shear_weights(source, weights_file, weight_type)
#         w_g = source['w']
#         print ('wg', w_g[:10])

#         return source['ra'] , source['dec'], source['e1'], source['e2'], R_g, w_g
    
    
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
        theta, correlation functions, and extra information
        """
        
        # minimum and maximum angular separation
        theta_min, theta_max = self.par.theta_lims
        
        treecorr_config = {'nbins': self.par.ang_nbins,
                         'min_sep': theta_min,
                         'max_sep': theta_max,
                         'sep_units': sep_units,
                         'bin_slop': self.par.bin_slop
                        }
        
        e1, e2, Rg, wg = params

        print('Num lenses = {:d}'.format(len(ra_l)))
        print('Num sources = {:d}'.format(len(ra_s)))
        print('Rg = {:f}'.format(Rg))
        print('wg = {}'.format(wg))
        print('Average e1 = {:f}'.format(np.average(e1, weights=wg)))
        print('Average e2 = {:f}'.format(np.average(e2, weights=wg)))
        
        # generate lens and source catalogs to correlate and process them
        ng = treecorr.NGCorrelation(treecorr_config)
        cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l, ra_units=units, dec_units=units, w=weights)
        cat_s = treecorr.Catalog(ra=ra_s, dec=dec_s, ra_units=units, dec_units=units,
                                 g1=(e1-np.average(e1, weights=wg)), 
                                 g2=(e2-np.average(e2, weights=wg)), 
                                 w=wg)
        ng.process(cat_l, cat_s, low_mem=low_mem)
        print('Done NG')

        # get theta, gamma_t
        theta = np.exp(ng.logr)
        gamma_t = ng.xi/Rg
        gammat_tot = np.copy(gamma_t)

        # get gamma_x
        gamma_x = ng.xi_im/Rg
        gammax_tot = np.copy(gamma_x)

        # generate randoms catalogs to correlate and process them
        if use_randoms or use_boosts:
            rg = treecorr.NGCorrelation(treecorr_config)
            cat_r = treecorr.Catalog(ra=ra_rand, dec=dec_rand, ra_units=units, dec_units=units)
            rg.process(cat_r, cat_s, low_mem=low_mem)
            print('Done RG')

        # get boost factors for gamma_t
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

        # random-point subtraction for gamma_t and gamma_x
        if use_randoms:
            gammat_rand = rg.xi/Rg
            gammax_rand = rg.xi_im/Rg
            xi_im_rand = rg.xi_im
            xi_npairs_rand = rg.npairs
            xi_weight_rand = rg.weight
        else:
            gammat_rand = np.zeros(len(theta))
            gammax_rand = np.zeros(len(theta))
            xi_im_rand = np.zeros(len(theta))
            xi_npairs_rand = np.zeros(len(theta))
            xi_weight_rand = np.zeros(len(theta))
        gammat_tot -= gammat_rand
        gammax_tot -= gammax_rand
        
        return (theta, 
                gammat_tot, gamma_t, gammat_rand, 
                gammax_tot, gamma_x, gammax_rand, 
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
        theta_min, theta_max = self.par.theta_lims
        
        treecorr_config = {'nbins': self.par.ang_nbins,
                         'min_sep': theta_min,
                         'max_sep': theta_max,
                         'sep_units': sep_units,
                         'bin_slop': self.par.bin_slop
                        }

        # ellipticity (e1,e2), R_gamma and weights
        e1, e2, Rg, wg = params

        print('Num lenses = {:d}'.format(len(ra_l)))
        print('Num sources = {:d}'.format(len(ra_s)))
        print('Rg = {:f}'.format(Rg))
        print('wg = {}'.format(wg))
        print('Average e1 = {:f}'.format(np.average(e1, weights=wg)))
        print('Average e2 = {:f}'.format(np.average(e2, weights=wg)))

        # generate lens and source catalogs to correlate and process them
        ng = treecorr.NGCorrelation(treecorr_config, var_method='jackknife')
        
        # generate randoms catalog to save jk_centers 
        cat_r = treecorr.Catalog(ra=ra_rand, dec=dec_rand, ra_units=units, dec_units=units, npatch=self.par.n_jck)
        cat_r.write_patch_centers('jk_centers')
        
        cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l, ra_units=units, dec_units=units, w=weights, patch_centers='jk_centers')
            
        cat_s = treecorr.Catalog(ra=ra_s, dec=dec_s, ra_units=units, dec_units=units,
                                 g1=(e1-np.average(e1, weights=wg)), 
                                 g2=(e2-np.average(e2, weights=wg)), 
                                 w=wg, patch_centers='jk_centers')
        
        ng.process(cat_l, cat_s, low_mem=low_mem)
        print('Done NG')

        # get theta, gamma_t
        theta = np.exp(ng.logr)
        gamma_t = ng.xi/Rg
        gammat_tot = np.copy(gamma_t)

        # get gamma_x
        gamma_x = ng.xi_im/Rg
        gammax_tot = np.copy(gamma_x)

        # generate randoms catalogs to correlate and process them
        if use_randoms or use_boosts:
            rg = treecorr.NGCorrelation(treecorr_config, var_method='jackknife')
            rg.process(cat_r, cat_s, low_mem=low_mem)
            print('Done RG')

        # get boost factors for gamma_t
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

        # random-point subtraction for gamma_t and gamma_x
        if use_randoms:
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

        # update correlations with responses to use in Jackknife mode
        ng.Rg = Rg*np.ones(len(theta))
            
        # generate fake treecorr correlation objects for lenses and randoms 
        # that hold the weights for the boost factor covariance calculations
        if use_boosts:
            # initialize NN correlations for single point cross lenses and randoms
            nn_lp = treecorr.NNCorrelation(nbins=self.par.ang_nbins, min_sep=1.e-3, max_sep=1.e5, 
                                           sep_units='arcmin', bin_slop=self.par.bin_slop, var_method='jackknife')
            nn_rp = treecorr.NNCorrelation(nbins=self.par.ang_nbins, min_sep=1.e-3, max_sep=1.e5, 
                                           sep_units='arcmin', bin_slop=self.par.bin_slop, var_method='jackknife')
            # catalog containing single point
            cat_p = treecorr.Catalog(ra=np.array([np.mean(ra_l)]), dec=np.array([np.mean(dec_l)]), ra_units=units, dec_units=units)
            # process catalogs and get objects 
            # containing the weights split into jackknife patches
            nn_lp.process(cat_l, cat_p, low_mem=low_mem)
            nn_rp.process(cat_r, cat_p, low_mem=low_mem)
            
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
        cov_jk_gt = treecorr.estimate_multi_cov(corrs, 'jackknife', func)

        # get gammax covariance
        if use_randoms:
            func = lambda corrs: corrs[0].xi_im/corrs[0].Rg - corrs[1].xi_im/corrs[1].Rg
            corrs = [ng,rg]
        else:
            func = lambda corrs: corrs[0].xi_im/corrs[0].Rg
            corrs = [ng]
        cov_jk_gx = treecorr.estimate_multi_cov(corrs, 'jackknife', func)

        # get boost factor covariance
        if use_boosts:
            func = lambda corrs: (corrs[0].weight/np.sum(corrs[2].weight)) / (corrs[1].weight/np.sum(corrs[3].weight))
            corrs = [ng,rg,nn_lp,nn_rp]
            cov_jk_boost = treecorr.estimate_multi_cov(corrs, 'jackknife', func)
        else:
            cov_jk_boost = np.zeros((len(theta),len(theta)))

        # get covariance of randoms points
        corrs = [rg]
        func = lambda corrs: corrs[0].xi/corrs[0].Rg
        cov_jk_rand = treecorr.estimate_multi_cov(corrs, 'jackknife', func)
        
        print('Done Covariance')

        return (theta, 
                gammat_tot, gamma_t, gammat_rand, 
                gammax_tot, gamma_x, gammax_rand, 
                cov_jk_gt, cov_jk_gx, cov_jk_boost, cov_jk_rand,
                ng.varxi,
                ng.xi_im, xi_im_rand, ng.npairs, xi_npairs_rand, ng.weight, xi_weight_rand, 
                Rg, sum_w_l, sum_w_r, boost)
    