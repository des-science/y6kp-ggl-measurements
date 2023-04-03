"""
----------------
Helper functions
----------------
"""
import numpy as np
import os
import sys
import yaml
import catalog_utils
import treecorr
from astropy.io import fits
from destest import destest

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

    def load_lens_Y3_maglim(self, data_file, zl_lims=None):
        """
        Loads lens galaxy data from file

        Options:
        - Y3 MagLim
        """
        # read data for lenses
        hdul = fits.open(data_file)
        ra = hdul[1].data['RA']
        dec = hdul[1].data['DEC']
        z = hdul[1].data['Z']
        w = hdul[1].data['W']

        # mask in redshift
        goodm = np.where( (z>zl_lims[0]) * (z<zl_lims[1]) )[0]
        
        return ra[goodm], dec[goodm], w[goodm]
    
    def load_randoms_Y3_maglim(self, data_file, zl_lims=None):
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
        
        return ra[goodm], dec[goodm]
    
    def load_source_metacal_5sels(self, data_file, zs_bin=None):
        """
        Loads source galaxy data

        Options:
        - Y3 Metacal
        """
        # read yaml file that defines all the catalog selections used
        params = yaml.load(open(data_file), Loader=yaml.FullLoader)
        params['param_file'] = data_file.split('/')[-1]

        source_selector, source_calibrator = catalog_utils.load_catalog(params, 'mcal', 'mcal', 
                                                                        params['source_group'], params['source_table'], 
                                                                        params['source_path'], 
                                                                        return_calibrator=destest.MetaCalib)

        gold_selector = catalog_utils.load_catalog(params, 'gold', 'mcal', params['gold_group'], params['gold_table'],
                                                    params['gold_path'], inherit=source_selector)

        pz_selector = catalog_utils.load_catalog(params, 'pz', 'mcal', params['pz_group'], 
                                                    params['pz_table'], params['pz_path'], inherit=source_selector)

        # dictionary with the unsheared version of each quantity with the 
        # selections from: unsheared, 1p, 1m, 2p, 2m
        source_5sels = {}
        source_5sels['unsheared'] = {}
        source_5sels['unsheared']['ra'] = [gold_selector.get_col('ra', uncut=True)[0][gold_selector.get_mask()[i]] for i in range(5)]
        source_5sels['unsheared']['dec'] = [gold_selector.get_col('dec', uncut=True)[0][gold_selector.get_mask()[i]] for i in range(5)]
        source_5sels['unsheared']['e1'] = [source_selector.get_col('e_1', uncut=True)[0][source_selector.get_mask()[i]] for i in range(5)]
        source_5sels['unsheared']['e2'] = [source_selector.get_col('e_2', uncut=True)[0][source_selector.get_mask()[i]] for i in range(5)]

        # dictionary with the 5 selections (1p, 1m, 2p, 2m), in which the 
        # quantities are obtained from sheared images or using fluxes measured
        # on sheared images
        source_5sels['sheared'] = {}
        source_5sels['sheared']['e1'] = source_selector.get_col('e_1')
        source_5sels['sheared']['e2'] = source_selector.get_col('e_2')
        # formerly z_g
        source_5sels['sheared']['som_bin'] = pz_selector.get_col('bhat')

        # dictionary with the unsheared version and selection only
        source = {}
        source['ra'] = source_5sels['unsheared']['ra'][0]
        source['dec'] = source_5sels['unsheared']['dec'][0]
        source['e1'] = source_5sels['unsheared']['e1'][0]
        source['e2'] = source_5sels['unsheared']['e2'][0]
        source['som_bin'] = source_5sels['sheared']['som_bin'][0]

        # masks from redshifts
        photoz_masks = [
                        (source_5sels['sheared']['som_bin'][i] >= zs_bin) & \
                        (source_5sels['sheared']['som_bin'][i] < zs_bin+1) for i in range(5)
                       ]
        maskz = photoz_masks[0]

        ra_s = source['ra'][maskz]
        dec_s = source['dec'][maskz]
        e1_s = source['e1'][maskz]
        e2_s = source['e2'][maskz]

        # calibration given the photoz bin
        R1,_,w_g = source_calibrator.calibrate('e_1', mask=photoz_masks)
        R2,_,w_g = source_calibrator.calibrate('e_2', mask=photoz_masks)
        R_g = 0.5*(R1+R2)

        return ra_s, dec_s, e1_s, e2_s, R_g, w_g
    
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
        theta_min, theta_max = self.par.theta_lims

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
                                 g1=(e1-np.average(e1, weights=wg)), 
                                 g2=(e2-np.average(e2, weights=wg)), 
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
        theta_min, theta_max = self.par.theta_lims

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
        if os.path.isfile('jk_centers'):
            cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l, ra_units=units, dec_units=units, w=weights, patch_centers='jk_centers')
        else:
            cat_r = treecorr.Catalog(ra=ra_rand, dec=dec_rand, ra_units=units, dec_units=units, npatch=self.par.n_jck)
            cat_r.write_patch_centers('jk_centers')
            cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l, ra_units=units, dec_units=units, w=weights, patch_centers='jk_centers')
        cat_s = treecorr.Catalog(ra=ra_s, dec=dec_s, ra_units=units, dec_units=units,
                                 g1=(e1-np.average(e1, weights=wg)), 
                                 g2=(e2-np.average(e2, weights=wg)), 
                                 w=wg, patch_centers='jk_centers')
        ng.process(cat_l, cat_s, low_mem=low_mem)

        # get theta, gammat
        theta = np.exp(ng.logr)
        gamma_t = ng.xi/Rg
        gammat_tot = gamma_t

        # get imaginary part of xi and gamma_x
        gamma_x = ng.xi_im/Rg
        gammax_tot = gamma_x

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
            nn_lp = treecorr.NNCorrelation(nbins=self.par.ang_nbins, min_sep=1.e-3, max_sep=1.e5, sep_units='arcmin', bin_slop=self.par.bin_slop, var_method='jackknife')
            nn_rp = treecorr.NNCorrelation(nbins=self.par.ang_nbins, min_sep=1.e-3, max_sep=1.e5, sep_units='arcmin', bin_slop=self.par.bin_slop, var_method='jackknife')
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
                func = lambda corrs: ( (corrs[0].weight/np.sum(corrs[2].weight)) / (corrs[1].weight/np.sum(corrs[3].weight)) 
                                        * corrs[0].xi/corrs[0].Rg 
                                        - corrs[1].xi/corrs[1].Rg )
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
        cov_jk_rand = treecorr.estimate_multi_cov(corrs, 'jackknife', func)

        return (theta, gamma_t, gammat_tot, gammat_rand, gamma_x, gammax_tot, gammax_rand, 
                cov_jk_gt, ng.varxi, cov_jk_boost, cov_jk_gx, cov_jk_rand,
                ng.xi_im, xi_im_rand, ng.npairs, xi_npairs_rand, ng.weight, xi_weight_rand, 
                Rg, sum_w_l, sum_w_r, boost)
    