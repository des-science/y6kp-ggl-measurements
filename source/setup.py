"""
----------------
Helper functions
----------------
"""
import numpy as np
import os
import yaml
import catalog_utils
import healpy as hp

class GGL_setup(object):

    def __init__(self):
        pass

    def load_lens(self, data_file):
        """
        Loads lens galaxy data from file

        Options:
        - Y3 MagLim
        - Y3 redMaGiC
        """
        # read data for lenses
        _data = np.loadtxt(data_file)
        ra  = _data[:,0] 
        dec = _data[:,1]
        z   = _data[:,2]
        
        return ra, dec, z
    
    def load_source_metacal_5sels(self, data_file):
        """
        Loads source galaxy data

        Options:
        - Y3 Metacal
        """
        # read yaml file that defines all the catalog selections used
        params = yaml.load(open(data_file), Loader=yaml.FullLoader)
        params['param_file'] = self.par.file_yaml

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

        return source, source_5sels, source_calibrator
    
    def mask_lens(self, ra_l, dec_l, z_l, zl_lims=None, mask_file=None, NSIDE=None, nest=False):
        """
        Define masks to apply to lens data sets
        """
        if mask_file is not None:
            # read mask
            mask = np.load(mask_file)

            theta = (90.0-dec_l)*np.pi/180.
            phi   = ra_l*np.pi/180.
            pix   = hp.ang2pix(NSIDE, theta, phi, nest=nest)
            goodm = np.where( (mask[pix]==1) * (z_l>zl_lims[0]) * (z_l<zl_lims[1]) )[0]
        else:
            goodm = np.where( (z_l>zl_lims[0]) * (z_l<zl_lims[1]) )[0]

        return ra_l[goodm], dec_l[goodm]
    
    def mask_source_metacal_5sels(self, ra_s, dec_s, e1_s, e2_s, source_5sels, calibrator, ra_jk=None, dec_jk=None, zs_bin=None):
        """
        Define masks to apply to source data sets
        """
        # masks from redshifts
        photoz_masks = [
                        (source_5sels['sheared']['som_bin'][i] >= zs_bin) & \
                        (source_5sels['sheared']['som_bin'][i] < zs_bin+1) for i in range(5)
                       ]
        maskz = photoz_masks[0]

        # apply just the redshift mask
        ra_s = ra_s[maskz]
        dec_s = dec_s[maskz]
        e1_s = e1_s[maskz]
        e2_s = e2_s[maskz]

        # calibration given the photoz bin
        R1,_,w_g = calibrator.calibrate('e_1', mask=photoz_masks)
        R2,_,w_g = calibrator.calibrate('e_2', mask=photoz_masks)
        R_g = 0.5*(R1+R2)

        # apply additional mask for Jackknife patch
        if (ra_jk is not None) and (dec_jk is not None):
            maskzjk = self.get_nearby_mask(ra_s, dec_s, ra_jk, dec_jk,  nside_mask=self.par.nside_nearby)
            ra_s = ra_s[maskzjk]
            dec_s = dec_s[maskzjk]
            e1_s = e1_s[maskzjk]
            e2_s = e2_s[maskzjk]
            w_g = w_g[maskzjk]

        return ra_s, dec_s, e1_s, e2_s, R_g, w_g
    