import psutil
import os
import numpy as np
import astropy.io.fits as pf
import pathos.multiprocessing as mp
from multiprocessing import Manager
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
import scipy.stats
import healpy as hp
import treecorr
import signal
import twopoint
from scipy import interpolate
import functions
import sys
import yaml
#sys.path.append('../../destest/')
#import destest
from destest import destest
import destest_functions
import ipdb
import h5py as h


def make_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


class GGL(object):
    """
    Basic class that has all the functions shared for several tests.
    """

    def __init__(self, basic, config, paths):
        self.basic = basic
        self.config = config
        self.paths = paths

    def get_nz(self,z):
        zbins = np.linspace(0,2.5,500)	
        zbinsc = zbins[:-1] + (zbins[1]-zbins[0])/2.
        nz, _ = np.histogram(z,zbins)
        nz = nz/float(nz.sum())
        return zbins, nz 

    def get_nz_weighted(self,z,w):
        zbins = np.linspace(0,2.5,500)	
        zbinsc = zbins[:-1] + (zbins[1]-zbins[0])/2.
        nz, _ = np.histogram(z,zbins, weights=w)
        nz = nz/float(nz.sum())
        return zbins, nz

    def get_responses(self,calibrator, mask=None, ind_responses=False):
        '''
        Function to deal with responses and make sure all things are good.
        It is a function so we can use it for several samples, for instance,
        for the non-tomographic sample, each of the tomographic bins,
        and other selections you could come up with.
        Makes several tests to ensure internal consistency and finally
        returns the relevant quantities.
        '''
        # These responses below do no include the weights,
        # they are just the full response (Rg+Rs) for each object, component and this selection
        # c is just zero for Y3, would be the additive response
        R11, c1, w1 = calibrator.calibrate('e_1', return_full=True, mask=mask) 
        R22, c2, w2 = calibrator.calibrate('e_2', return_full=True, mask=mask)
        R = (R11+R22)/2.

        # since (if) the weights are the same we can just keep one column
        assert ((w1==w2).all())
        w = w1

        # Get the total mean response directly
        R11_mean = calibrator.calibrate('e_1', mask=mask)[0] # this is a scalar
        R22_mean = calibrator.calibrate('e_2', mask=mask)[0]
        R_mean = np.mean([R11_mean, R22_mean])

        # And check everything holds
        assert(R11_mean == np.average(R11,weights=w))
        assert(R22_mean == np.average(R22,weights=w))

        # Lets get the quantity we need to weight our redshift
        # distributions with: weights*fullresponse(Rgamma+Rs)
        wR = calibrator.calibrate('e_1', return_wRgS=True, mask=mask)
        wRg = calibrator.calibrate('e_1', return_wRg=True, mask=mask)

        # Responses are tricky! So again lets check it all makes sense
        assert(np.isclose(wR,w*R).all())
        Rg = wRg/w
        Rg_mean = np.average(Rg, weights=w)
        Rs_mean = R_mean-Rg_mean
        assert (np.isclose(Rs_mean,(np.average(wR/w, weights=w) - Rg_mean)))

        # Let's make some prints
        print('From xcorr: R = %0.4f'%R_mean) 
        print('From xcorr: Rg = %0.4f'%Rg_mean)
        print('From xcorr: Rs = %0.4f'%Rs_mean)
        print('From xcorr: R11 = %0.4f'%R11_mean) 
        print('From xcorr: R22 = %0.4f'%R22_mean)

        # Save in two separate dictionaries, one to save to file, one needed later in the code
        save = {}
        save['R_mean'] = R_mean
        save['Rg_mean'] = Rg_mean
        save['Rs_mean'] = Rs_mean
        save['R11_mean'] = R11_mean
        save['R22_mean'] = R22_mean

        res = {}
        res['R_mean'] = R_mean
        if ind_responses:
            # include the individual responses to perform the scale dependence averaging test
            res['R'] = R 
        #res['wR'] = wR # would be needed for constructing the N(z) 
        res['w'] = w #weights but for convienence we load them here
        
        return res, save

    @profile
    def load_metacal(self, reduce_mem=False):
        """
        Loads metacal data for Y3 catalog using h5 interface.
        Combines with Gold and BPZ.
        Returns: two dictionaries for the sources, with all the relevant columns and the calibrator
                 in order to be used for other functions too. 
                 These are the two dictionaries it returns, which are used for different things:

        - source_5sels: It is a nested dictionary with 'sheared' and 'unsheared' quantities. 
                        *'unsheared': Dictionary with the unsheared version of each 
                                      quantity with the selections from: unsheared, 1p, 1m, 2p, 2m.
                                      The unsheared version of a quantity means that has been measured on
                                      images which are unsheared, or using fluxes obtained from unsheared images.
                                      This dictionary is useful to compute the selection response manually.

                        *'sheared': Dictionary with the 5 selections (1p, 1m, 2p, 2m), in which the 
                                    quantities are obtained from sheared images or using fluxes measured
                                    on sheared images. 
                                    This dictionary is useful to make further selections besides the baseline selection,
                                    for instance to select of redshift bins, size or S/N, which you need to do also 
                                    using the sheard quanties to compute the selection response. 
        
        - source: Simplest one. Dictionary with the baseline unsheared selection and quantites.
        - reduce_mem: If true, loads only essential columns to run the measurement, to save ram memory.
        
        """

        # Read yaml file that defines all the catalog selections used
        params = yaml.load(open(self.paths['yaml']))
        params['param_file'] = self.paths['yaml']

        # Source catalog
        source_selector, source_calibrator = destest_functions.load_catalog(
            params, 'mcal', 'mcal', params['source_group'], params['source_table'], params['source_path'],
            return_calibrator=destest.MetaCalib)

        # Gold catalog
        gold_selector = destest_functions.load_catalog(
            params, 'gold', 'mcal', params['gold_group'], params['gold_table'], params['gold_path'],
            inherit=source_selector)

        # BPZ (or DNF) catalog, depending on paths in cats.yaml file (exchange bpz and dnf)
        #bpz_selector = destest_functions.load_catalog(
        #    params, 'bpz', 'mcal', params['bpz_group'], params['bpz_table'], params['bpz_path'], inherit=source_selector)

        # SOM PZ to split in bins:
        som_selector = destest_functions.load_catalog(
            params, 'pz', 'mcal', params['pz_group'], params['pz_table'], params['pz_path'], inherit=source_selector) 
        
        # Dictionary with the unsheared version of each quantity with the selections from: unsheared, 1p, 1m, 2p, 2m. 
        source_5sels = {}
        if not reduce_mem:
            source_5sels['unsheared'] = {}
            source_5sels['unsheared']['ra'] = [gold_selector.get_col('ra', uncut=True)[0][gold_selector.get_mask()[i]] for i
                                               in range(5)]
            source_5sels['unsheared']['dec'] = [gold_selector.get_col('dec', uncut=True)[0][gold_selector.get_mask()[i]] for
                                                i in range(5)]
            source_5sels['unsheared']['e1'] = [source_selector.get_col('e_1', uncut=True)[0][source_selector.get_mask()[i]]
                                               for i in range(5)]
            source_5sels['unsheared']['e2'] = [source_selector.get_col('e_2', uncut=True)[0][source_selector.get_mask()[i]]
                                               for i in range(5)]

        # Dictionary with the 5 selections (1p, 1m, 2p, 2m), in which the quantities are obtained from sheared images or using fluxes measured on sheared images.   
        source_5sels['sheared'] = {}
        source_5sels['sheared']['som_bin'] = som_selector.get_col('bhat')
        if not reduce_mem:
            source_5sels['sheared']['e1'] = source_selector.get_col('e_1')
            source_5sels['sheared']['e2'] = source_selector.get_col('e_2')
            source_5sels['sheared']['snr'] = source_selector.get_col('snr')
            source_5sels['sheared']['size'] = source_selector.get_col('T')
            source_5sels['sheared']['bpz_mean'] = bpz_selector.get_col('zmean_sof')
            source_5sels['sheared']['bpz_zmc'] = bpz_selector.get_col('zmc_sof')

        # Dictionary with the unsheared version and selection only:
        source = {}
        if reduce_mem:
            source['ra'] = gold_selector.get_col('ra')[0]
            source['dec'] = gold_selector.get_col('dec')[0]
            source['e1'] = source_selector.get_col('e_1')[0]
            source['e2'] = source_selector.get_col('e_2')[0]
            
        if not reduce_mem:
            source['ra'] = source_5sels['unsheared']['ra'][0]
            source['dec'] = source_5sels['unsheared']['dec'][0]
            source['e1'] = source_5sels['unsheared']['e1'][0]
            source['e2'] = source_5sels['unsheared']['e2'][0]
            source['psf_e1'] = source_selector.get_col('psf_e1')[0]
            source['psf_e2'] = source_selector.get_col('psf_e2')[0]
            source['snr'] = source_5sels['sheared']['snr'][0]
            source['size'] = source_5sels['sheared']['size'][0]
            source['bpz_mean'] = source_5sels['sheared']['bpz_mean'][0]
            source['bpz_zmc'] = source_5sels['sheared']['bpz_zmc'][0]
            source['som_bin'] = source_5sels['sheared']['som_bin'][0]

        # lets deal with responses, check function above
        dict_responses, dict_to_file = self.get_responses(source_calibrator)

        # Lets include the responses dictionary into the source dictionary
        source = dict(source, **dict_responses)

        # Lets save some of the mean responses values
        # might be useful to put the numbers in the papers        
        import json
        with open(self.get_path_test_allzbins() + 'responses_untomographic', 'w') as file:
            file.write(json.dumps(dict_to_file))
        # can be loaded with: json.load(open('test'))
         
        return source, source_5sels, source_calibrator

    @profile
    def load_metacal_bin(self, source, source_5sels, calibrator, bin_low, CalCa, reduce_mem = False, ind_responses=False):
        """
        source: dictionary containing relevant columns for the sources, with the baseline selection applied already.
        source_5sels: dictionary containing relevant columns for the sources, with the baseline selection applied already,
                     for each of the 5 selections 1p, 1m, 2p, 2m. 
        calibrator: class to compute the response. Taken from baseline selection.
        bin_low, bin_high: limits to select the tomographic bin. Can be 0,1,2,3.
        Obtains 5 masks (unsheared, sheared 1p, 1m, 2p, 2m) to obtain the new selection response.
        Returns: Source dictionary masked with the unsheared mask and with the mean response updated.
        """

        photoz_masks = [
            (source_5sels['sheared']['som_bin'][i] >= bin_low ) & (source_5sels['sheared']['som_bin'][i] <= bin_high)
            for i in range(5)]
        source_bin = {}
        source_bin['ra'] = source['ra'][photoz_masks[0]]
        source_bin['dec'] = source['dec'][photoz_masks[0]]
        source_bin['e1'] = source['e1'][photoz_masks[0]]
        source_bin['e2'] = source['e2'][photoz_masks[0]]
        if not reduce_mem:
            source_bin['psf_e1'] = source['psf_e1'][photoz_masks[0]]
            source_bin['psf_e2'] = source['psf_e2'][photoz_masks[0]]
            source_bin['snr'] = source['snr'][photoz_masks[0]]
            source_bin['size'] = source['size'][photoz_masks[0]]
            source_bin['bpz_mean'] = source['bpz_mean'][photoz_masks[0]]
            source_bin['bpz_zmc'] = source['bpz_zmc'][photoz_masks[0]]

        # lets deal with responses, check function above
        dict_responses, dict_to_file = self.get_responses(calibrator, mask=photoz_masks, ind_responses=ind_responses)

        # Lets include the responses dictionary into the source dictionary
        source_bin = dict(source_bin, **dict_responses)

        # Lets save some of the mean responses values
        # might be useful to put the numbers in the papers        
        import json
        with open(self.get_path_test_allzbins() + 'responses_%d_%d'%(bin_low, bin_high), 'w') as file:
            file.write(json.dumps(dict_to_file))
        # can be loaded with: json.load(open('test'))
        return source_bin

    def subtract_mean_shear(self, cat):
        """
        Subtracts the weighted mean shear from a source catalog, 
        from each component.
        """
        mean_e1 = np.average(cat['e1'], weights=cat['w'])
        mean_e2 = np.average(cat['e2'], weights=cat['w'])

        cat['e1'] -= mean_e1
        cat['e2'] -= mean_e2

        print('Subtracting mean shear.')
        print('Mean e1:', mean_e1)
        print('Mean e2:', mean_e2)

        return cat, np.array([mean_e1, mean_e2])
     
    def get_lens(self, lens):
        """
        Given a lens sample, returns ra, dec, jk and weight, in case it exists.
        """

        ra_l = lens['ra']
        dec_l = lens['dec']
        jk_l = lens['jk']
        try:
            w_l = lens['w']
            print('Weights found in foreground catalog.')
        except:
            print('There are no identified weights for the foreground sample.')
            #w_l = np.ones(len(ra_l))
            w_l = None


        if not self.config['lens_w']:
            print('Running in mode with no LSS weights. They are set to one.')
            w_l = np.ones(len(ra_l)) 

        return ra_l, dec_l, jk_l, w_l

    def get_source(self, source):
        """
        Given a source sample, returns ra, dec, and weight, in case it exists.
        """

        ra_s = source['ra']
        dec_s = source['dec']
        try:
            w = source['w']
            print('Weights found in source catalog.')
        except:
            print('There are no identified weights for the sources.')
            w = np.ones(len(ra_s))

        return ra_s, dec_s, w

    @profile
    def run_treecorr_jackknife(self, lens, source, type_corr):
        """
        Function that runs treecorr for a given lens and source sample,
        and a given configuration, and paralalizes the measurements
        in different jackknife patches for the lenses and randoms.
        Returns the measurements for each jackknife patch.
        type_corr: string, type of correlation, i.e. NG, NN, NK_Rgamma
        NG for gammat, NN for wtheta, NK for scalar quantities, like
        the responses. The string after NK_ is used to load the column.
        """
        assert type_corr == 'NG' or type_corr == 'NN' or 'NK' in type_corr, 'This type_corr of correlation is not accepted by this function.'
        parent_id = os.getpid()

        def worker_init():
            """
            Function to easily kill the run when there is multiprocessing.
            """

            def sig_int(signal_num, frame):
                print('signal: %s' % signal_num)
                parent = psutil.Process(parent_id)
                for child in parent.children():
                    if child.pid != os.getpid():
                        print("killing child: %s" % child.pid)
                        child.kill()
                print("killing parent: %s" % parent_id)
                parent.kill()
                print("suicide: %s" % os.getpid())
                psutil.Process(os.getpid()).kill()

            signal.signal(signal.SIGINT, sig_int)
        @profile
        def run_jki_NG(jk):
            """
            Function we use for mutiprocessing.
            jk: Region we run in each core.
            For NG correlations.
            """
            print(jk)
            ra_l_jk = ra_l[jk_l==jk]
            dec_l_jk = dec_l[jk_l==jk]
            if w_l is not None:
                print('There are weights inside NG')
                w_l_jk = w_l[jk_l==jk]

            if jk == 0: print('Doing NG correlation.')
            corr = treecorr.NGCorrelation(nbins=self.config['nthbins'], min_sep=self.config['thlims'][0],
                                          max_sep=self.config['thlims'][1], sep_units='arcmin',
                                          bin_slop=self.config['bslop'], num_threads=self.basic['num_threads'])

            if len(ra_l_jk) > 1:

                if len(ra_l_jk)%2 == 0:  # If len array is even, make it odd before computing the median.
                    pixjk = hp.ang2pix(nside, ((90.0 - np.median(np.append(dec_l_jk, 0)))*np.pi/180.),
                                       np.median(np.append(ra_l_jk, 0))*np.pi/180.0)
                if len(ra_l_jk)%2 == 1:  # If len array is odd, keep it odd.
                    pixjk = hp.ang2pix(nside, ((90.0 - np.median(dec_l_jk))*np.pi/180.),
                                       np.median(ra_l_jk)*np.pi/180.)

                pixsjk = hp.get_all_neighbours(nside, pixjk)
                pixsjk = np.append(pixsjk, pixjk)
                bool_s = np.in1d(pix, pixsjk)
                #del pixsjk
                #del pixjk

                if w_l is not None:
                    print('a) Using weights')
                    cat_l = treecorr.Catalog(ra=ra_l_jk, dec=dec_l_jk, w=w_l_jk, ra_units='deg', dec_units='deg')
                else:
                    print('b) Not using weights')
                    cat_l = treecorr.Catalog(ra=ra_l_jk, dec=dec_l_jk, ra_units='deg', dec_units='deg')
                
                if self.config['source_only_close_to_lens']:
                    cat_s = treecorr.Catalog(ra=ra_s[bool_s], dec=dec_s[bool_s], g1=e1[bool_s], g2=e2[bool_s],
                                             w=w[bool_s], ra_units='deg', dec_units='deg')
                else:
                    cat_s = treecorr.Catalog(ra=ra_s, dec=dec_s, g1=e1, g2=e2, w=w, ra_units='deg', dec_units='deg')
                            
                corr.process(cat_l, cat_s)

                if jk == 0: theta.append(np.exp(corr.logr))
                gts[jk].append(corr.xi)
                gxs[jk].append(corr.xi_im)
                errs[jk].append(np.sqrt(np.abs(corr.varxi)))
                weights[jk].append(corr.weight)
                npairs[jk].append(corr.npairs)

            else:
                if jk == 0: theta.append(np.exp(corr.logr))
                zeros = np.zeros(self.config['nthbins'])
                gts[jk].append(zeros)
                gxs[jk].append(zeros)
                errs[jk].append(zeros)
                weights[jk].append(zeros)
                npairs[jk].append(zeros)

        def run_jki_NK(jk):
            """
            Function we use for mutiprocessing.
            jk: Region we run in each core.
            For NK correlations.
            """

            ra_l_jk = ra_l[jk_l == jk]
            dec_l_jk = dec_l[jk_l == jk]
            w_l_jk = w_l[jk_l == jk]

            if jk == 0: print('Doing NK correlation with variable %s.' % type_corr[3:])
            corr = treecorr.NKCorrelation(nbins=self.config['nthbins'], min_sep=self.config['thlims'][0],
                                          max_sep=self.config['thlims'][1], sep_units='arcmin',
                                          bin_slop=self.config['bslop'], num_threads=self.basic['num_threads'])

            if len(ra_l_jk) > 1:

                if len(ra_l_jk) % 2 == 0:  # If len array is even, make it odd before computing the median.
                    pixjk = hp.ang2pix(nside, ((90.0 - np.median(np.append(dec_l_jk, 0))) * np.pi / 180.),
                                       np.median(np.append(ra_l_jk, 0)) * np.pi / 180.0)
                if len(ra_l_jk) % 2 == 1:  # If len array is odd, keep it odd.
                    pixjk = hp.ang2pix(nside, ((90.0 - np.median(dec_l_jk)) * np.pi / 180.),
                                       np.median(ra_l_jk) * np.pi / 180.)

                pixsjk = hp.get_all_neighbours(nside, pixjk)
                pixsjk = np.append(pixsjk, pixjk)
                bool_s = np.in1d(pix, pixsjk)
 
                cat_l = treecorr.Catalog(ra=ra_l_jk, dec=dec_l_jk, w=w_l_jk, ra_units='deg', dec_units='deg')
                
                if self.config['source_only_close_to_lens']:
                    cat_s = treecorr.Catalog(ra=ra_s[bool_s], dec=dec_s[bool_s], k=scalar[bool_s], w=w[bool_s],
                                                 ra_units='deg', dec_units='deg')
                else:
                    cat_s = treecorr.Catalog(ra=ra_s, dec=dec_s, k=scalar, w=w, ra_units='deg', dec_units='deg')

                corr.process(cat_l, cat_s)

                if jk == 0: theta.append(np.exp(corr.logr))
                xi_nks[jk].append(corr.xi)
                weights[jk].append(corr.weight)
                npairs[jk].append(corr.npairs)

            else:
                if jk == 0: theta.append(np.exp(corr.logr))
                zeros = np.zeros(self.config['nthbins'])
                xi_nks[jk].append(corr.xi)
                weights[jk].append(zeros)
                npairs[jk].append(zeros)

        ra_l, dec_l, jk_l, w_l = self.get_lens(lens)
        if not self.config['lens_w']:
            print('Checking weights of lens sample are one:', w_l)
        ra_s, dec_s, w = self.get_source(source)

        # We choose a healpix grid to select sources around each lens patch, to reduce memory
        # Choose nside for this grid here
        nside = self.config['nside']
        theta = (90.0 - dec_s) * np.pi / 180.
        phi = ra_s * np.pi / 180.
        pix = hp.ang2pix(nside, theta, phi)

        if type_corr == 'NG':
            e1 = source['e1']
            e2 = source['e2']
        if 'NK' in type_corr:
            print(type_corr, type_corr[3:])
            scalar = source['%s' % type_corr[3:]]

        manager = Manager()
        theta = manager.list()
        weights = [manager.list() for x in range(self.config['njk'])]
        npairs = [manager.list() for x in range(self.config['njk'])]
        if type_corr == 'NG':
            gts = [manager.list() for x in range(self.config['njk'])]
            gxs = [manager.list() for x in range(self.config['njk'])]
            errs = [manager.list() for x in range(self.config['njk'])]
            
        if 'NK' in type_corr:
            xi_nks = [manager.list() for x in range(self.config['njk'])]

        if self.basic['pool']: 
            p = mp.Pool(self.basic['Ncores'], worker_init)
            if type_corr == 'NG':
                p.map(run_jki_NG, range(self.config['njk']))
            if 'NK' in type_corr:
                p.map(run_jki_NK, range(self.config['njk']))
            p.close()

        if not self.basic['pool']:
            for jk in range(self.config['njk']):
                if type_corr == 'NG':
                    run_jki_NG(jk)
                if 'NK' in type_corr:
                    run_jki_NK(jk)

        def reshape_manager(obj):
            return (np.array(list(obj))).reshape(self.config['njk'], self.config['nthbins'])

        print('Reshaping manager... theta')
        theta = (np.array(list(theta))).reshape(self.config['nthbins'])
        print('Reshaping manager... weights')
        weights = reshape_manager(weights)
        print('Reshaping manager... npairs')
        npairs = reshape_manager(npairs)

        if type_corr == 'NG':
            print('Reshaping manager... gts, gx, err')
            gts = reshape_manager(gts)
            gxs = reshape_manager(gxs)
            errs = reshape_manager(errs)
            return theta, gts, gxs, errs, weights, npairs

        if 'NK' in type_corr:
            print('returning NK')
            xi_nks = reshape_manager(xi_nks)
            return theta, xi_nks, weights, npairs

    def run_nk_no_jackknife(self, lens, source, scalar):
        """
        Uses TreeCorr to compute the NK correlation between lens and source.
        Used to compute scale dependece responses.
        scalar: array with scalar values of some quantity of which we want to compute the
               average in angular bins.
        Returns theta and R_nk.
        """

        ra_l, dec_l, jk_l, w_l = self.get_lens(lens)
        ra_s, dec_s, w = self.get_source(source)

        nk = treecorr.NKCorrelation(nbins=self.config['nthbins'], min_sep=self.config['thlims'][0],
                                    max_sep=self.config['thlims'][1], sep_units='arcmin', bin_slop=self.config['bslop'], num_threads=self.config['num_threads'])

        cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l, w=w_l, ra_units='deg', dec_units='deg')
        cat_s = treecorr.Catalog(ra=ra_s, dec=dec_s, k=scalar, w=w, ra_units='deg', dec_units='deg')
        nk.process(cat_l, cat_s)
        theta = np.exp(nk.logr)
        R_nk = nk.xi

        return theta, R_nk

    def save_runs(self, path_test, theta, gts, gxs, errs, weights, npairs, random_bool):
        """
        Function to save the measurements in each single jk patch.
        """
        rp = '_rp' if random_bool else ''
        np.savetxt(path_test + 'theta' + rp, theta, header='theta [arcmin]')
        np.savetxt(path_test + 'gts' + rp, gts, header='gts for all jackknife regions')
        np.savetxt(path_test + 'gxs' + rp, gxs, header='gxs for all jackknife regions')
        np.savetxt(path_test + 'errs' + rp, errs, header='errs for all jackknife regions')
        np.savetxt(path_test + 'weights' + rp, weights, header='weights for all jackknife regions')
        np.savetxt(path_test + 'npairs' + rp, npairs, header='npairs for all jackknife regions')

    def load_runs(self, path_test, random_bool):
        """
        Function to load the measurements in each single jk patch.
        """
        rp = '_rp' if random_bool else ''
        theta = np.loadtxt(path_test + 'theta' + rp)
        gts = np.loadtxt(path_test + 'gts' + rp)
        gxs = np.loadtxt(path_test + 'gxs' + rp)
        errs = np.loadtxt(path_test + 'errs' + rp)
        weights = np.loadtxt(path_test + 'weights' + rp)
        npairs = np.loadtxt(path_test + 'npairs' + rp)

        return theta, gts, gxs, errs, weights, npairs

    def compute_boost_factor_jackknife(self, jk_l, jk_r, wnum, wnum_r, w_l):
        """
        Computes the boost factor for a given set of weights for the lenses and randoms.
        Uses Eq. from Sheldon et al. (2004)
        It does so for each N-1 jk regions.
        jk_l: number of assigned jk region per each galaxy
        jk_r: number of assigned jk region per random point
        wnum: array of weights for the lenses, for all jk regions -1, each time different region.
        wnum_r: array of weights for the randoms, for all jk regions -1, each time different region. 
        w_l: weights assigned to each lens galaxy.
        """

        bf_jk = []
        jks = np.arange(self.config['njk'])
        ratio_n_jk = np.zeros(self.config['njk'])
        for jk in jks:
            ratio_n_jk[jk] = float((jk_r != jk).sum()) / ((jk_l != jk)*w_l).sum()
            bf_jk.append(ratio_n_jk[jk] * wnum[jk] / wnum_r[jk])

        bf_jk = np.array(bf_jk)
        return bf_jk


    def compute_boost_factor_exact(self, wnum, wnum_r, lens, random):
        """
        Computes the boost factor for a given set of weights for the lenses and randoms.
        Uses Eq. from Sheldon et al. (2004)
        wnum: array of weights as a function of scale for the lens-source correlations, as output of TreeCorr.
        wnum_r: similar array for random-source correlations.
        lens: lens catalog including the single object weights.
        random: random catalog which might include weights or not.
        """
        try:
            n_lens = np.sum(lens['w'])
        except:
            n_lens = len(lens['ra'])

        try:
            w_r = random['w']
            n_ran = np.sum(w_r)
        except:
            n_ran = len(random['ra'])

        bf = np.array((n_ran/n_lens)*(wnum/wnum_r))
        
        return bf

    def numerators_jackknife(self, gts, gxs, ws):
        """
        Inputs:
        - gts: Tangential shear as output from TreeCorr for each single JK patch.
        - gxs: Same object for the cross-component.
        - ws: Same object for the weights, as they come out of TreeCorr, as a function of scale.

        Returns:
        - The numerators of gt, gx and weights, as the sum of each of the N-1 patches, leaving
        one out at a time. These quantites are combined later to get the JK covariance matrix and 
        mean tangential shear from JK. The exact tangential shear (as if done without JK) is obtaned
        using the function numerators_exact. 

        Note: This function assumes the run was done with process, not process_cross treecorr functions.
        Assuming that means that gts are already normalized by the weights, but here we want to 
        'unnormalize' this quantity so it can be combined between different JK regions. 
        (Advantage of process vs process_cross is that process allows you to save the shape noise err)

        """
        # Obtain the numerators, see note above.
        gts = gts * ws
        gxs = gxs * ws

        # Construct array taking the sum of all gt from each jk patch except for one patch, different each time
        gt_num = np.array([np.sum(np.delete(np.array(gts), i, axis=0), axis=0) for i in range(len(gts))])
        gx_num = np.array([np.sum(np.delete(np.array(gxs), i, axis=0), axis=0) for i in range(len(gts))])
        w_num = np.array([np.sum(np.delete(np.array(ws), i, axis=0), axis=0) for i in range(len(gts))])

        return gt_num, gx_num, w_num

    def numerators_exact(self, gts, gxs, ws):
        """
        Inputs:
        - gts: Tangential shear as output from TreeCorr for each single JK patch.
        - gxs: Same object for the cross-component.
        - ws: Same object for the weights, as they come out of TreeCorr, as a function of scale.

        Returns:
        - The numerators of gt, gx and weights, as the sum the N patches. This function recovers the 
        tangential shear as if it was not done with JK. 

        Note: This function assumes the run was done with process, not process_cross treecorr functions.
        Assuming that means that gts are already normalized by the weights, but here we want to 
        'unnormalize' this quantity so it can be combined between different JK regions. 
        (Advantage of process vs process_cross is that process allows you to save the shape noise err)
        """
        # Obtain the numerators, see note above.
        gts = gts * ws
        gxs = gxs * ws

        # Sum result from each JK patch to combine the mesurements
        gt_num = np.sum(gts, axis=0) 
        gx_num = np.sum(gxs, axis=0) 
        w_num = np.sum(ws, axis=0) 

        return gt_num, gx_num, w_num

    def process_run(self, exact, all, theta, path_test, end):
        """
        From the jackknife measurements in all jackknife regions but all, constructs covariance, mean and stats.
        Saves them into file.
        exact: gt or gx as a sum of all JK regions, exact result as if we were not doing JK.
        all: gt_all or gx_all.
        theta: in arcmin.
        path_test: where to save the files.
        end: string to save the files: gt, gx, randoms, etc.
        """
        mean = np.mean(all, axis=0)
        cov = functions.covariance(all, mean)
        err = np.sqrt(np.diag(cov))
        chi2 = np.dot(mean.T, np.dot(np.linalg.inv(cov), mean))
        ndf = len(mean)

        std = np.sqrt((len(all) - 1.)) * np.std(all, axis=0)

        # Hartlap factor
        N = len(all)  # number of jackknife regions
        p = len(mean)  # number of angular bins
        factor = (N - p - 2) / float(N - 1)
        chi2_hartlap = chi2 * factor

        stats = np.array([chi2_hartlap, ndf])

        np.savetxt(path_test + 'mean_JK_%s' % end, zip(theta, mean, err), header='th, %s, err_%s' % (end, end))
        np.savetxt(path_test + 'cov_%s' % end, cov)
        np.savetxt(path_test + 'all_%s' % end, all,
                   header='%s (sum of %s from all patches except for one, different each time)' % (end, end))
        np.savetxt(path_test + 'null_chi2_%s' % end, stats.reshape(1, stats.shape[0]),
                   fmt='%0.1f  %d', header='chi2_hartlap  ndf')
        np.savetxt(path_test + '%s' % end, zip(theta, exact, err), header='th, %s, err_%s' % (end, end))

    def get_chi2(self, path_test, end):
        """
        Load chi2 and degrees of freedom for a certain lens and source bin and test.
        end: string like: gt, gx, randoms.
        """
        chi2, ndf = np.loadtxt(path_test + 'null_chi2_%s' % end, unpack=True)
        return chi2, ndf

    def save_plot(self, name_plot):
        """
        Saves plot in pdf and png and creates directory to save it if it doesn't exist.
        """

        make_directory(self.paths['plots_config'])
        plt.savefig(self.paths['plots_config'] + '%s.pdf' % name_plot, bbox_inches='tight')
        plt.savefig(self.paths['plots_config'] + '%s.png' % name_plot, bbox_inches='tight', dpi=400)
        print('Plot saved in:', self.paths['plots_config'] + '%s.pdf' % name_plot)

    def load_sims(self):
        """
	    Loads the simulation measurements and covariance used to fit an amplitude to the data measurements.
        Function used in TestSizeSNR and TestSysMaps.
	    """
        sims = np.loadtxt(self.paths['sims'], unpack=True, usecols=(1,))
        cov_sims = np.loadtxt(self.paths['cov_sims'])
        return sims, cov_sims

    def ratio_from_sims(self, theta, gtl_all, gth_all, sims, cov_sims):
        """
	    Computes the ratio of measurements of the splits (high/low) by computing the ratio of the corresponding amplitudes with respect to the simulation measurements.
        Function used in TestSizeSNR and TestSysMaps.
	    """
        ratioA_all = np.zeros(len(gth_all))
        chi2fit_h = np.zeros(len(gth_all))
        chi2fit_l = np.zeros(len(gth_all))

        mask = (theta >= self.source_nofz_pars['thetamin'])
        for i in range(len(gth_all)):
            Ah, chi2fit_h[i], _ = functions.minimize_chi2_fit_amplitude((cov_sims[mask].T)[mask].T, gth_all[i][mask],
                                                                        sims[mask])
            Al, chi2fit_l[i], _ = functions.minimize_chi2_fit_amplitude((cov_sims[mask].T)[mask].T, gtl_all[i][mask],
                                                                        sims[mask])
            ratioA_all[i] = Ah / Al

        ratioA_mean = np.mean(ratioA_all, axis=0)
        COV_A = (len(ratioA_all) - 1) * np.cov(ratioA_all, bias=True)
        err_A = np.sqrt(COV_A)
        return ratioA_mean, err_A

    def save_2pointfile(self, string):
        """
        Save the correlation function measurements (i.e. gammat, boost factors, etc),
        N(z)'s and jackknife covariance into the 2point format file.
        Creates a:
        - SpectrumMeasurement object: In which the 2PCF measurements are saved. 
        - Kernel object: The N(z)'s are here. 
        - CovarianceMatrixInfo object: In which the jackknife covariance is saved.

        Then, it builds the TwoPointFile objet from the above objects,
        and saves it to a file.
        """
             
        # adapted from save_2pt, removed unused options
        def get_scales(x_min, x_max, nbins):
            """
            Get scales
            """
            log_lims = np.linspace(np.log(x_min), np.log(x_max), nbins+1 )
            lims = np.exp(log_lims)
            xmin = lims[:-1]
            xmax = lims[1:]
            mids = (2./3.) * (xmax**3 - xmin**3) / (xmax**2 - xmin**2)
            return lims, mids
            
        # Preparing spectrum
        length = self.config['nthbins'] * len(self.zbins['lbins']) * len(self.zbins['sbins'])
        values = np.zeros(length, dtype=float)
        bin1 = np.zeros(length, dtype=int)
        bin2 = np.zeros_like(bin1)
        angular_bin = np.zeros_like(bin1)
        angle = np.zeros_like(values)
        angle_min = np.zeros_like(values)
        angle_max = np.zeros_like(values)
        dv_start = 0
        cov = np.zeros((length, length))
        # getting the edges of the angular bins
        # Also replacing theta to be exactly the same as theory lines
        edges, mids = get_scales(self.config['thlims'][0], self.config['thlims'][1], self.config['nthbins'])
        
        for l in range(0, len(self.zbins['lbins'])):
            for s in range(len(self.zbins['sbins'])):
                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                theta, xi, xi_err = np.loadtxt(path_test + '%s' % string, unpack=True)
                cov_ls = np.loadtxt(path_test + 'cov_%s' % string)

                bin_pair_inds = np.arange(dv_start, dv_start + self.config['nthbins'])
                values[bin_pair_inds] = xi
                bin1[bin_pair_inds] = l + 1
                bin2[bin_pair_inds] = s + 1
                angular_bin[bin_pair_inds] = np.arange(self.config['nthbins'])
                angle[bin_pair_inds] = mids
                # now these needed in the twopoint files for blinding
                angle_min[bin_pair_inds] = edges[0:-1]
                angle_max[bin_pair_inds] = edges[1:]
                dv_start += self.config['nthbins']
                cov[bin_pair_inds[0]:bin_pair_inds[-1] + 1, bin_pair_inds] = cov_ls

        # Preparing N(z) 
        if 'data' in self.basic['mode']:
            # if it is a TwoPointFile use the code below
            file_lens_nz = twopoint.TwoPointFile.from_fits(self.paths['lens_nz'])
            lens_nz = file_lens_nz.get_kernel('nz_lens')
            file_source_nz = twopoint.TwoPointFile.from_fits(self.paths['source_nz'])
            source_nz = file_source_nz.get_kernel('nz_source')

            print('Saving TwoPointFile with lens N(z) from %s'%(self.paths['lens_nz']))
            print('Saving TwoPointFile with source N(z) from %s'%(self.paths['source_nz']))

            '''
            # If instead of a TwoPointFile we only have the number densities objects
            # use this code below.
            lensfile = pf.open(self.paths['lens_nz'])
            lens_nz = twopoint.NumberDensity.from_fits(lensfile[2])
            assert (lens_nz.name=='nz_lens')
            
            sourcefile = pf.open(self.paths['source_nz'])
            source_nz = twopoint.NumberDensity.from_fits(lensfile[1])
            assert (source_nz.name=='nz_source')
            '''
            
        if 'gt' in string:
            gammat = twopoint.SpectrumMeasurement('gammat', (bin1, bin2),
                                                  (twopoint.Types.galaxy_position_real,
                                                   twopoint.Types.galaxy_shear_plus_real),
                                                  ('nz_lens', 'nz_source'), 'SAMPLE', angular_bin, values,
                                                  angle=angle, angle_unit='arcmin',
                                                  angle_min=angle_min, angle_max=angle_max) 

            cov_mat_info = twopoint.CovarianceMatrixInfo('COVMAT', ['gammat'], [length], cov)
            gammat_twopoint = twopoint.TwoPointFile([gammat], [lens_nz, source_nz], windows=None, covmat_info=cov_mat_info)
            twopointfile_unblind = self.get_twopointfile_name(string)

            # Remove file if it exists already because to_fits function doesn't overwrite
            if os.path.isfile(twopointfile_unblind):
                os.system('rm %s' % (twopointfile_unblind))

            gammat_twopoint.to_fits(twopointfile_unblind)

        if string == 'boost_factor':
            # saves the boost factors independently in a fits file
            boost_factor = twopoint.SpectrumMeasurement('boost_factor', (bin1, bin2),
                                                        (twopoint.Types.galaxy_position_real,
                                                         twopoint.Types.galaxy_shear_plus_real),
                                                        ('nz_lens', 'nz_source'), 'SAMPLE', angular_bin, values,
                                                        angle=angle, angle_unit='arcmin')

            cov_mat_info = twopoint.CovarianceMatrixInfo('COVMAT', ['boost_factor'], [length], cov)

            print('Saving TwoPointFile')
            boost_factor_twopoint = twopoint.TwoPointFile([boost_factor], [lens_nz, source_nz], windows=None, covmat_info=cov_mat_info)
            save_path = os.path.join(self.get_path_test_allzbins() + '%s_twopointfile.fits'%string)

            # Remove file if it exists already because to_fits function doesn't overwrite                                                                                                                                                                                    
            if os.path.isfile(save_path):
                os.system('rm %s' % (save_path))

            boost_factor_twopoint.to_fits(save_path)

        if string == 'randoms':
            # saves the random points measurements in a fits file
            random_points = twopoint.SpectrumMeasurement('random_points', (bin1, bin2),
                                                        (twopoint.Types.galaxy_position_real,
                                                         twopoint.Types.galaxy_shear_plus_real),
                                                        ('nz_lens', 'nz_source'), 'SAMPLE', angular_bin, values,
                                                        angle=angle, angle_unit='arcmin')

            cov_mat_info = twopoint.CovarianceMatrixInfo('COVMAT', ['random_points'], [length], cov)

            print('Saving TwoPointFile')
            random_points_twopoint = twopoint.TwoPointFile([random_points], [lens_nz, source_nz], windows=None, covmat_info=cov_mat_info)
            save_path = os.path.join(self.get_path_test_allzbins() + '%s_twopointfile.fits'%string)

            # Remove file if it exists already because to_fits function doesn't overwrite                                                                                                                                                                                    
            if os.path.isfile(save_path):
                os.system('rm %s' % (save_path))

            random_points_twopoint.to_fits(save_path)

        if string == 'gx':
            # saves the random points measurements in a fits file
            random_points = twopoint.SpectrumMeasurement('gx', (bin1, bin2),
                                                        (twopoint.Types.galaxy_position_real,
                                                         twopoint.Types.galaxy_shear_plus_real),
                                                        ('nz_lens', 'nz_source'), 'SAMPLE', angular_bin, values,
                                                        angle=angle, angle_unit='arcmin')

            cov_mat_info = twopoint.CovarianceMatrixInfo('COVMAT', ['gx'], [length], cov)

            print('Saving TwoPointFile')
            random_points_twopoint = twopoint.TwoPointFile([random_points], [lens_nz, source_nz], windows=None, covmat_info=cov_mat_info)
            save_path = os.path.join(self.get_path_test_allzbins() + '%s_twopointfile.fits'%string)

            # Remove file if it exists already because to_fits function doesn't overwrite                                                                                                                                                                                    
            if os.path.isfile(save_path):
                os.system('rm %s' % (save_path))

            random_points_twopoint.to_fits(save_path)

            
            
    def load_data_or_sims(self):
        '''
        Loads and returns lens, randoms and sources, either from data or from simulations.
        '''

        if self.basic['mode'] == 'data':
            # load Y3 sources here
            source_all, source_all_5sels, calibrator = self.load_metacal(reduce_mem=True)
            return source_all, source_all_5sels, calibrator

        if self.basic['mode'] == 'mice':
            lens_all = pf.getdata(self.paths['lens_mice'])
            random_all = pf.getdata(self.paths['randoms_mice'])
            return lens_all, random_all

        if self.basic['mode'] == 'buzzard':
            if 'redmagic' in self.config['lens_v']:
                lens_all = pf.getdata(self.paths['lens_buzzard'])
                random_all = pf.getdata(self.paths['randoms_buzzard'])
                source_all = self.load_buzzard()
                return lens_all, random_all, source_all
            
            if 'maglim' in self.config['lens_v']:
                source_all = self.load_buzzard()
                return source_all



class Measurement(GGL):
    """
    Subclass that runs the gglensing measurement for all the lens-source bin combinations.
    Includes:
    - Mean response calculation.
    - Random points subtraction.
    - Jackknife covariance calculation.
    - Boost factors calculation.
    """

    def __init__(self, basic, config, paths, zbins, plotting):
        GGL.__init__(self, basic, config, paths)
        self.zbins = zbins
        self.plotting = plotting

    def get_path_test(self, lbin, sbin):
        return os.path.join(self.paths['runs_config'], 'measurement', lbin + '_' + sbin) + '/'

    def get_path_test_allzbins(self):
        return os.path.join(self.paths['runs_config'], 'measurement') + '/'

    def get_twopointfile_name(self, string):
        return os.path.join(self.get_path_test_allzbins() + '%s_twopointfile.fits' % string)

    def get_y1_mask(self, ra, dec):
        hdul = pf.open('/global/project/projectdirs/des/ggl/lens_cats/Y1/' + 'DES_Y1A1_3x2pt_redMaGiC_MASK_HPIX4096RING.fits')
        y1 = hdul[1].data
        fracgood = y1['FRACGOOD']
        hpix = y1['HPIX']
        map = np.zeros(hp.nside2npix(4096))+hp.UNSEEN
        map[hpix] = fracgood
        # convert ra, dec to pix
        pix_d = hp.ang2pix(4096, ra, dec, lonlat=True,nest=0)
        mask = map[pix_d]>0.8
        return mask
        
    def sel_lens_zbin(self, lens_all, lbin):
        """
        Select lens bins depending on whether its data or sims.
        Input: 
        - lens_all: dictionary with all the relevant columns of the lenses: ra, dec, z, w, jk.
        - lbin: which lens bin to select, e.g., 'l1'. Redshift edges defined in info.py file.
        Returns:
        - Lens: Selected reshift bin object with the same info.

        """
        # Decide how to bin the lenses
        if 'ztrue' in self.config['zllim_v']:
            print('Using ztrue to bin the lenses.')
            lens = lens_all[(lens_all['ztrue'] > self.zbins[lbin][0]) & (lens_all['ztrue'] < self.zbins[lbin][1])]
        else:
            print('Using z to bin the lenses.')
            lens = lens_all[(lens_all['z'] > self.zbins[lbin][0]) & (lens_all['z'] < self.zbins[lbin][1])]

        print('Length lens', lbin, len(lens['ra']), self.zbins[lbin][0], self.zbins[lbin][1])

        if self.basic['mode'] == 'buzzard':
            zbins, nz_l = self.get_nz(lens['ztrue'])		    
            np.savetxt(self.get_path_test_allzbins()+'/nzs/'+'nz_%s'%lbin,nz_l)

        if self.basic['mode'] == 'mice':
            zbins, nz_l = self.get_nz(lens['z'])		    
            np.savetxt(self.get_path_test_allzbins()+'/nzs/'+'nz_%s'%lbin,nz_l)

        return lens
 

    def sel_random_zbin(self, random_all, lbin):
        """
        Similar function as sel_lens_zbin for random points.
        """

        if 'data' in self.basic['mode']:
            random = random_all[(random_all['z'] > self.zbins[lbin][0]) & (random_all['z'] < self.zbins[lbin][1])]
        if self.basic['mode']  == 'buzzard' or self.basic['mode'] == 'mice':
            random = random_all[(random_all['z'] > self.zbins[lbin][0]) & (random_all['z'] < self.zbins[lbin][1])]
        if self.basic['mode']  == 'mice_old':
            random = random_all[l*len(random_all)/len(self.zbins['lbins']):(l+1)*len(random_all)/len(self.zbins['lbins'])]
        return random

    @profile
    def run(self):
        make_directory(self.get_path_test_allzbins()+'/nzs/')

        if self.basic['mode'] == 'data':
            mean_shears = []

        if self.basic['mode'] == 'mice':
            lens_all, random_all = self.load_data_or_sims()
            
        if self.basic['mode'] == 'buzzard':
            print('Running on Buzzard..')
            if 'redmagic' in self.config['lens_v']:
                lens_all, random_all, source_all = self.load_data_or_sims()
            if 'maglim' in self.config['lens_v']:
                source_all = self.load_data_or_sims()

        for sbin in self.zbins['sbins']:

            print('Running measurement for source %s.' % sbin)
    
            if self.basic['mode'] == 'data':
                source_all, source_all_5sels, calibrator = self.load_data_or_sims()
                print('Source bin:', self.zbins[sbin][0], self.zbins[sbin][1])
                source = self.load_metacal_bin(source_all, source_all_5sels, calibrator, bin_low=self.zbins[sbin][0], bin_high=self.zbins[sbin][1], reduce_mem=True, ind_responses=True)

                print('Before..')
                print(np.average(source['R'],weights=source['w']))
                print(source['R_mean'])
                
                if self.basic['area'] == 'y1': 
                    y1mask = self.get_y1_mask(source['ra'], source['dec'])
                    for key in source:
                        if key is not 'R_mean':
                            source[key] = source[key][y1mask]

                source['R_mean'] = np.average(source['R'],weights=source['w'])
                                                    
                R_mean = source['R_mean']
                print('Length source', sbin, len(source['ra']))
                print(R_mean)

                source, mean_shear = self.subtract_mean_shear(source)
                print(mean_shear)
                mean_shears.append(mean_shear)

                # Delete these huge dictionaries to save memory while running treecorr (otherwise crashes, giving bus error)
                del source_all
                del source_all_5sels
                del calibrator
                del source['R']

            if self.basic['mode'] == 'data_y1sources':
                source = pf.getdata(self.paths['y1_sources'] + 'metacal_sel_sa%s.fits'%sbin[1])
                print('Length source', sbin, len(source['ra']))
                R_s_bins = [0.0072, 0.014, 0.0098, 0.014] # from Y1 gglensing paper
                print("Sel response:", R_s_bins[int(sbin[1])-1])
                print("Rgamma:", np.average(source['Rgamma']))
                R_mean = np.average(source['Rgamma']) + R_s_bins[int(sbin[1])-1]
                print("Total response:", R_mean)

            elif self.basic['mode'] == 'mice_y1sources':
                """
                In this case there are no responses, so we set it to one.
                """
                R_mean = 1.
                #source = source_all[(source_all['z'] > self.zbins[sbin][0]) & (source_all['z'] < self.zbins[sbin][1])]
                hdu = pf.open(self.paths['mice_y1sources'] + 'mice2_shear_fullsample_downsampled025_bins%s.fits'%sbin[-1])
                mice = hdu[1].data

                if np.isnan(mice['ra']).any():
                    print('There is a nan in the source catalog, in ra.')
                    mask_nan = np.isnan(mice['ra'])==False
                    print('There are %d nans.'%(len(mice['ra']) - mask_nan.sum()))
                    
                if np.isnan(mice['dec']).any():
                    print('There is a nan in the source catalog, in dec.')

                source = {}
                source['ra'] = mice['ra']
                source['dec'] = mice['dec']
                source['e1'] = -mice['e1'] # flip e1 sign
                source['e2'] = mice['e2']
                source['ztrue'] = mice['ztrue']
                source['w'] = mice['w']
                
                print('Length source', sbin, len(source['ra']))
                print('np.std(e1)', np.std(source['e1']))
                print('np.std(e2)', np.std(source['e2']))
                zbins, nz_s = self.get_nz_weighted(source['ztrue'], source['w'])		    
                np.savetxt(self.get_path_test_allzbins()+'/nzs/'+'zbins',zbins,header='zbin limits')
                np.savetxt(self.get_path_test_allzbins()+'/nzs/'+'nz_%s'%sbin,nz_s)

            elif self.basic['mode'] == 'mice':
                """
                In this case there are no responses, so we set it to one.
                """
                R_mean = 1.
                hdu = pf.open(self.paths['mice_y1sources'] + 'new_dnf_sof_mice_Feb3rd_zbin%s_imax23.fits'%sbin[-1])
                mice = hdu[1].data

                if np.isnan(mice['ra_gal']).any():
                    print('There is a nan in the source catalog, in ra.')
                    mask_nan = np.isnan(mice['ra_gal'])==False
                    print('There are %d nans.'%(len(mice['ra_gal']) - mask_nan.sum()))
                    
                if np.isnan(mice['dec_gal']).any():
                    print('There is a nan in the source catalog, in dec.')

                def add_intrinsic_ellipticity(g1, g2, e1, e2):
                    # Shear field                                                                                                                                   
                    g = g1 + 1j*g2
                    # Intrinsic ellipticity                                                                                                                         
                    e = e1 + 1j*e2
                    res = (e+g)/(1.+np.conjugate(g)*e)
                    return res.real, res.imag

                e1,e2 = add_intrinsic_ellipticity(mice['gamma1'], mice['gamma2'], mice['eps1'], mice['eps2'])

                if np.isnan(e1).any():
                    print('There is a nan in the source catalog, in e1.')

                mask_nan = np.isnan(e1)==False

                if np.isnan(e1).any():
                    print('There are %d nans.'%(len(e1) - mask_nan.sum()))

                source = {}
                source['ra'] = mice['ra_gal'][mask_nan]
                source['dec'] = mice['dec_gal'][mask_nan]
                source['e1'] = -e1[mask_nan] # flip e1 sign
                source['e2'] = e2[mask_nan]
                #source['ztrue'] = mice['z_cgal'][mask_nan]


                print('Length source', sbin, len(source['ra']))
                print('np.std(e1)', np.std(source['e1']))
                print('np.std(e2)', np.std(source['e2']))
                print('np.std(ra)', np.std(source['ra']))
                print('np.std(dec)', np.std(source['dec']))
                zbins, nz_s = self.get_nz(mice['z_cgal'][mask_nan])		    
                np.savetxt(self.get_path_test_allzbins()+'/nzs/'+'zbins',zbins,header='zbin limits')
                np.savetxt(self.get_path_test_allzbins()+'/nzs/'+'nz_%s'%sbin,nz_s)

                del mice
                del hdu
                    
                    
            elif self.basic['mode'] == 'buzzard':
                """
                In this case there are no responses, so we set it to one.
                """
                R_mean = 1.
                source = {}
                for k in source_all.keys():
                    source[k] = source_all[k][(source_all['zbin'] >= self.zbins[sbin][0]) & (source_all['zbin'] <= self.zbins[sbin][1])]
                    print('Length source', sbin, len(source['ra']))
                    print('np.std(e1)', np.std(source['e1']))
                    print('np.std(e2)', np.std(source['e2']))
                    
                zbins, nz_s = self.get_nz(source['ztrue'])		    
                np.savetxt(self.get_path_test_allzbins()+'/nzs/'+'zbins',zbins,header='zbin limits')
                np.savetxt(self.get_path_test_allzbins()+'/nzs/'+'nz_%s'%sbin,nz_s)

                for l, lbin in enumerate(self.zbins['lbins']):
                    path_test = self.get_path_test(lbin, sbin)

                    if os.path.exists(path_test + 'gt_boosted'):
                        print('Measurements for this bin already exist. SKIPPING!')

                    else:
                        print('Running measurement for lens %s.' % lbin)
                        make_directory(path_test)

                    # Lenses run
                    if 'data' in self.basic['mode']:
                        lens_all = pf.getdata(self.paths['lens'])

                    if self.basic['mode'] == 'buzzard' and 'maglim' in self.config['lens_v']:
                        lens = pf.getdata(self.paths['buzzard'] + 'lens_%d.fits'%(l+1))
                        print('Loading lenses from:', self.paths['buzzard'] + 'lens_%d.fits'%(l+1))
                    else:
                        lens = self.sel_lens_zbin(lens_all, lbin)
                    print('Number of lenses in bin %s'%lbin, len(lens['ra']))
                        
                    if 'data' in self.basic['mode']:
                        # remove this object to save memory before calling Treecorr.
                        del lens_all
                    if np.isnan(lens['ra']).any():
                        print('There is a nan in the lens catalog, in ra.')
                    if np.isnan(lens['dec']).any():
                        print('There is a nan in the lens catalog, in dec.')

                    theta, gts, gxs, errs, weights, npairs = self.run_treecorr_jackknife(lens, source, 'NG')
                    self.save_runs(path_test, theta, gts, gxs, errs, weights, npairs, random_bool=False)
                    #theta, gts, gxs, errs, weights, npairs = self.load_runs(path_test, random_bool=False)
                    gtnum_jk, gxnum_jk, wnum_jk = self.numerators_jackknife(gts, gxs, weights)
                    gtnum_ex, gxnum_ex, wnum_ex = self.numerators_exact(gts, gxs, weights)
                    make_directory(self.get_path_test_allzbins()+'/weights/')
                    np.savetxt(self.get_path_test_allzbins()+'/weights/'+'w_%s_%s'%(lbin, sbin),wnum_ex,header='weights (sum of all JK regions)')

                    # Randoms run
                    if 'data' in self.basic['mode']:
                        # Randoms run
                        random_all = pf.getdata(self.paths['randoms'])
                    if self.basic['mode'] == 'buzzard' and 'maglim' in self.config['lens_v']:
                        random = pf.getdata(self.paths['buzzard'] + 'random_%d.fits'%(l+1))
                        print('Loading randoms from:', self.paths['buzzard'] + 'random_%d.fits'%(l+1))
                    else:
                        random = self.sel_random_zbin(random_all, lbin)
                    print('Number of randoms in bin %s'%lbin, len(random['ra']))

                    if 'data' in self.basic['mode']:
                        # remove this object to save memory before calling Treecorr.
                        del random_all

                    theta, gts, gxs, errs, weights, npairs = self.run_treecorr_jackknife(random, source, 'NG')
                    self.save_runs(path_test, theta, gts, gxs, errs, weights, npairs, random_bool=True)
                    #theta, gts, gxs, errs, weights, npairs = self.load_runs(path_test, random_bool=True)
                    gtnum_jk_r, gxnum_jk_r, wnum_jk_r = self.numerators_jackknife(gts, gxs, weights)
                    gtnum_ex_r, gxnum_ex_r, wnum_ex_r = self.numerators_exact(gts, gxs, weights)

                    # Combine final estimator for JK covariance and mean, with and without boost factors
                    gt_all = (gtnum_jk / wnum_jk) / R_mean - (gtnum_jk_r / wnum_jk_r) / R_mean
                    gx_all = (gxnum_jk / wnum_jk) / R_mean - (gxnum_jk_r / wnum_jk_r) / R_mean

                    try:
                        w_l = lens['w']
                        print('Weights found in foreground catalog.')
                    except:
                        print('There are no identified weights for the foreground sample.')
                        w_l = np.ones(len(lens['ra']))

                    bf_all = self.compute_boost_factor_jackknife(lens['jk'], random['jk'], wnum_jk, wnum_jk_r, w_l)
                    gt_all_boosted = bf_all*(gtnum_jk / wnum_jk)/R_mean - (gtnum_jk_r / wnum_jk_r) / R_mean 

                    # Combine final estimator for exact results (as done without JK)
                    gt = (gtnum_ex / wnum_ex) / R_mean - (gtnum_ex_r / wnum_ex_r) / R_mean
                    gx = (gxnum_ex / wnum_ex) / R_mean - (gxnum_ex_r / wnum_ex_r) / R_mean                        
                    bf = self.compute_boost_factor_exact(wnum_ex, wnum_ex_r, lens, random)
                    gt_boosted = bf*(gtnum_ex/wnum_ex)/R_mean - (gtnum_ex_r/wnum_ex_r)/R_mean 

                    self.process_run(gt, gt_all, theta, path_test, 'gt')
                    self.process_run(gx, gx_all, theta, path_test, 'gx')
                    self.process_run((gtnum_ex_r/wnum_ex_r)/R_mean, (gtnum_jk_r/wnum_jk_r)/R_mean, theta, path_test, 'randoms')
                    self.process_run(bf, bf_all, theta, path_test, 'boost_factor')
                    self.process_run(gt_boosted, gt_all_boosted, theta, path_test, 'gt_boosted')

            elif self.basic['mode'] == 'data':
                np.savetxt(self.get_path_test_allzbins()+'mean_shears', mean_shears, header = 'Mean_e1 Mean_e2 (for each redshift bin)')

            else:
                print('mode: "%s" is not implemented'%self.basic['mode'])

    def save_spectrum_measurement_file(self):
        """
        Save the boost factors into a SpectrumMeasurment file, from the twopoint code.
        Caveat: cannot save the errors, use TwoPointFile instead (see save_2pointfile function above).
        """

        bf_length = self.config['nthbins'] * len(self.zbins['lbins']) * len(self.zbins['sbins'])
        bf_values = np.zeros(bf_length, dtype=float)
        bin1 = np.zeros(bf_length, dtype=int)
        bin2 = np.zeros_like(bin1)
        angular_bin = np.zeros_like(bin1)
        angle = np.zeros_like(bf_values)
        dv_start = 0

        for l in range(0, len(self.zbins['lbins'])):
            for s in range(len(self.zbins['sbins'])):
                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                theta, bf, bf_err = np.loadtxt(path_test + 'boost_factor', unpack=True)

                bin_pair_inds = np.arange(dv_start, dv_start + self.config['nthbins'])
                bf_values[bin_pair_inds] = bf
                bin1[bin_pair_inds] = l
                bin2[bin_pair_inds] = s
                angular_bin[bin_pair_inds] = np.arange(self.config['nthbins'])
                angle[bin_pair_inds] = theta
                dv_start += self.config['nthbins']

        boost_factors = twopoint.SpectrumMeasurement('boost_factors', (bin1, bin2),
                                                     (twopoint.Types.galaxy_position_real,
                                                      twopoint.Types.galaxy_shear_plus_real),
                                                     ['no_nz', 'no_nz'], 'SAMPLE', angular_bin, bf_values,
                                                     angle=angle, angle_unit='arcmin', extra_cols=None)

        path_save = self.get_path_test_allzbins()
        filename = 'boost_factors.fits'
        (boost_factors.to_fits()).writeto(path_save + filename)

    def plot_NO_BLINDED(self):
        """"
        Makes plot of the fiducial measurement for all redshift bins.
        """

        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        cmap = self.plotting['cmap']
        cmap_step = 0.25
        title_source = self.plotting['catname']

        # Figure
        fig, ax = plt.subplots(2, 3, figsize=(10, 6), sharey=True, sharex=False, gridspec_kw={'height_ratios': [1, 1]})
        fig.subplots_adjust(hspace=0.0, wspace=0.00)

        for l in range(0, len(self.zbins['lbins'])):

            # To iterate between the three columns and two lines
            j = 0 if l < 3 else 1
            ax[j][l % 3].axvspan(self.config['thlims'][0] * 0.8, self.plotting['th_limit'][l], color='gray', alpha=0.2)
            for s in range(len(self.zbins['sbins'])):

                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                if os.path.isfile(path_test + 'gt'):
                    th, gt, err = np.loadtxt(path_test + 'gt', unpack=True)

                    mask_neg = gt < 0
                    mask_pos = gt > 0

                    # chi2, ndf = self.get_chi2(path_test, 'gt')
                    ax[j][l % 3].errorbar(th[mask_neg] * (1 + 0.05 * s), -gt[mask_neg], err[mask_neg], fmt='.',
                                          mfc='None',
                                          mec=plt.get_cmap(cmap)(cmap_step * s),
                                          ecolor=plt.get_cmap(cmap)(cmap_step * s), capsize=2)
                    ax[j][l % 3].errorbar(th[mask_pos] * (1 + 0.05 * s), gt[mask_pos], err[mask_pos], fmt='.',
                                          color=plt.get_cmap(cmap)(cmap_step * s),
                                          mec=plt.get_cmap(cmap)(cmap_step * s), label=self.plotting['redshift_s'][s],
                                          capsize=2)

                    ax[j][l % 3].set_xlim(self.config['thlims'][0] * 0.8, self.config['thlims'][1] * 1.2)
                    ax[j][l % 3].set_xscale('log')
                    ax[j][l % 3].set_yscale('log')

                    ax[j][l % 3].text(0.5, 0.9, self.plotting['redshift_l'][l], horizontalalignment='center',
                                      verticalalignment='center', transform=ax[j][l % 3].transAxes, fontsize=12)
                    ax[j][l % 3].text(0.5, 0.93, self.plotting['titles_redmagic'][l], horizontalalignment='center',
                                      verticalalignment='center', transform=ax[j][l % 3].transAxes, fontsize=12)


                    #if l % 3 > 0:  # In case we want to keep labels on the left y-axis
                    #ax[j][l % 3].yaxis.set_ticklabels([])  # to remove the ticks labels
                    #if l < 2:
                        #ax[0][l].xaxis.set_ticklabels([])  # to remove the ticks labels


                    ax[j][l % 3].set_xlabel(r'$\theta$ [arcmin]', size='large')
                    ax[j][0].set_ylabel(r'$\gamma_t (\theta)$', size='large')
                    # ax[j][l % 3].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%d$'))

                    """
                    # Chi2
                    ax[j][l%3].text(0.25,0.3,r'Null $\chi^2$/ndf',
                                    horizontalalignment='center', verticalalignment='center', transform=ax[j][l%3].transAxes, fontsize = 10)
                    ax[j][l%3].text(0.25,0.23 -0.06*s,r'$%0.1f/%d$'%(chi2, ndf),
                             horizontalalignment='center', verticalalignment='center', transform=ax[j][l%3].transAxes, fontsize = 12, color = plt.get_cmap(cmap)(cmap_step*s))
                    """

        ax[1][0].set_ylim(10 ** (-6), 0.999 * 10 ** (-2))
        ax[1][1].set_ylim(10 ** (-6), 0.999 * 10 ** (-2))

        #handles, labels = ax[0][0].get_legend_handles_labels()
        fig.delaxes(ax[1, 2])
        #ax[1][1].legend(handles[::-1], labels[::-1], frameon=True, fancybox = True,prop={'size':12}, numpoints = 1, loc='center left', bbox_to_anchor=(1, 0.5))
        #ax[0][0].legend(frameon=False, fancybox=True, prop={'size': 12}, numpoints=1, loc='center',
        #                bbox_to_anchor=(2.45, -0.52))
        fig.suptitle(title_source, fontsize=16)
        fig.subplots_adjust(top=0.93)
        self.save_plot('plot_measurement')

    def load_twopointfile(self, string):
        '''
        Loads TwoPointFile and returns it.
        '''
        filename = self.get_twopointfile_name(string)
        if self.basic['blind']:
            gammat_file = twopoint.TwoPointFile.from_fits('%s_BLINDED.fits' % filename[:-5])
        else:
            gammat_file = twopoint.TwoPointFile.from_fits('%s.fits' % filename[:-5])
        return gammat_file

    def compute_sn_ratio(self, string):
        '''
        Compute S/N ratio using null chi2. 
        S/N = sqrt(null chi2 - ndof)
        Uses full jackknife covariance.
        '''

        gammat_file = self.load_twopointfile(string)
        gammat = gammat_file.spectra[0].value
        cov = gammat_file.covmat

        COV = np.array(cov)
        INVCOV = np.linalg.inv(COV)
        null_chi2 = np.mat(gammat) * INVCOV * np.transpose(np.mat(gammat))
        dof = gammat.shape[0]
        sn = np.sqrt(null_chi2 - dof)

        path_save = self.get_path_test_allzbins()
        print('S/N of the full measurements %s:'%string, float(sn))
        np.savetxt(path_save + 'sn_ratio_full_measurements_%s' % string, sn, fmt='%0.5g',
                   header='S/N ratio computed with JK covariance, S/N = sqrt(null chi2 - ndof)')

    def plot_from_twopointfile(self, string):

        """
        Makes plot of the fiducial measurement for all redshift bins, from a twopoint file, like Y1 style.
        It also uses the twopoint plotting functions to plot the measurements in a different style,
        the covariance and the N(z)'s.
        Useful to plot the blinded measurements (now the default). 
        """

        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        gammat_file = self.load_twopointfile(string)

        gammat = gammat_file.spectra[0]
        pairs = gammat.bin_pairs
        npairs = len(pairs)
        # It starts with 1, not with 0
        bins_l = np.transpose(pairs)[0]
        bins_s = np.transpose(pairs)[1]
        nbins_l = np.max(bins_l)
        nbins_s = np.max(bins_s)
        assert len(self.zbins['lbins']) == nbins_l, 'Number of lens bins in info does not match with the one in the two-point file.'
        assert len(self.zbins['sbins']) == nbins_s, 'Number of source bins in info does not match with the one in the two-point file.'


        if self.plotting['use_cmap']:
            cmap = self.plotting['cmap']
            cmap_step = 0.25
            colors = []
            for s in range(len(self.zbins['sbins'])):
                colors.append(plt.get_cmap(cmap)(cmap_step * s))
        else:
            colors = self.plotting['colors']
        title_source = self.plotting['catname']

        # Figure
        fig, ax = plt.subplots(2, 3, figsize=(10, 6), sharey=True, sharex=False, gridspec_kw={'height_ratios': [1, 1]})
        fig.subplots_adjust(hspace=0.0, wspace=0.00)

        for l in range(0, len(self.zbins['lbins'])):
            # To iterate between the three columns and two lines
            j = 0 if l < 3 else 1
            ax[j][l % 3].axvspan(self.config['thlims'][0] * 0.8, self.plotting['th_limit'][l], color='gray', alpha=0.2)

            for s in range(len(self.zbins['sbins'])):
                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                th, gt = gammat.get_pair(l+1, s+1)
                err = gammat.get_error(l+1, s+1)

                mask_neg = gt < 0
                mask_pos = gt > 0

                ax[j][l % 3].errorbar(th[mask_neg] * (1 + 0.05 * s), -gt[mask_neg], err[mask_neg], fmt='.', mfc='None',
                                        mec=colors[s], ecolor=colors[s], capsize=2)
                ax[j][l % 3].errorbar(th[mask_pos] * (1 + 0.05 * s), gt[mask_pos], err[mask_pos], fmt='.',
                                        color=colors[s],
                                        mec=colors[s], label=self.plotting['redshift_s'][s], capsize=2)

                ax[j][l % 3].set_xlim(self.config['thlims'][0]*0.8, self.config['thlims'][1]*1.2)
                ax[j][l % 3].set_xscale('log')
                ax[j][l % 3].set_yscale('log')
                if self.config['thlims'][0] ==2.5:
                    ax[j][l % 3].set_ylim(10**(-6), 0.999*10**(-2))
                if self.config['thlims'][0] ==0.1:
                    ax[j][l % 3].set_ylim(10**(-6), 10**(-1))
                ax[j][l % 3].text(0.5, 0.9, self.plotting['redshift_l'][l], horizontalalignment='center',
                                    verticalalignment='center', transform=ax[j][l % 3].transAxes, fontsize=12)
                #ax[j][l % 3].text(0.5, 0.93, self.plotting['titles_redmagic'][l], horizontalalignment='center',
                #                  verticalalignment='center', transform=ax[j][l % 3].transAxes, fontsize=12)

                #if l % 3 > 0:  # In case we want to keep labels on the left y-axis
                ax[j][l % 3].yaxis.set_ticklabels([])  # to remove the ticks labels
                if l < 2:
                    ax[0][l].xaxis.set_ticklabels([])  # to remove the ticks labels

                ax[j][l % 3].set_xlabel(r'$\theta$ [arcmin]', size='large')
                ax[j][0].set_ylabel(r'$\gamma_t (\theta)$', size='large')
                ax[j][l % 3].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%d$'))

                """
                # Chi2
                ax[j][l%3].text(0.25,0.3,r'Null $\chi^2$/ndf',
                                horizontalalignment='center', verticalalignment='center', transform=ax[j][l%3].transAxes, fontsize = 10)
                ax[j][l%3].text(0.25,0.23 -0.06*s,r'$%0.1f/%d$'%(chi2, ndf),
                            horizontalalignment='center', verticalalignment='center', transform=ax[j][l%3].transAxes, fontsize = 12, color = plt.get_cmap(cmap)(cmap_step*s))
                """

    
        # ax[1][1].legend(handles[::-1], labels[::-1], frameon=True, fancybox = True,prop={'size':12}, numpoints = 1, loc='center left', bbox_to_anchor=(1, 0.5))
        if self.config['lens_v'] == 'redmagic':
            fig.delaxes(ax[1, 2])
            ax[0][0].legend(frameon=False, fancybox=True, prop={'size': 12}, numpoints=1, loc='center',bbox_to_anchor=(2.45, -0.52))
        fig.suptitle(title_source, fontsize=16)
        fig.subplots_adjust(top=0.93)
        if self.basic['blind']:
            blind = '_BLINDED'
        else:
            blind = ''

        self.save_plot('plot_measurement%s_%s'%(blind, string))

        gammat_file.plots(self.paths['plots_config'] + '%s_twopointfile%s'%(string, blind),
                          blind_yaxis=self.basic['blind'], latex=self.plotting['latex'])

        # Use twopoint library to make the rest of the plots
    def plot_boostfactors(self):

        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        cmap = self.plotting['cmap']
        fig, ax = plt.subplots(len(self.zbins['sbins']), len(self.zbins['lbins']), figsize=(12.5, 9.375), sharey=True, sharex=True)
        fig.subplots_adjust(hspace=0.1, wspace=0.1)
        c1 = plt.get_cmap(cmap)(0.)

        for l in range(0, len(self.zbins['lbins'])):
            for s in range(len(self.zbins['sbins'])):
                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                theta, bf, bf_err = np.loadtxt(path_test + 'boost_factor', unpack=True)

                ax[s][l].axhline(y=1, ls=':', color='k', alpha=1)
                ax[s][l].errorbar(theta, bf, bf_err, fmt='.', color=c1, mec=c1, capsize=2)
                ax[s][l].set_xscale('log')
                ax[s][l].set_xlim(self.config['thlims'][0], self.config['thlims'][1])
                ax[s][l].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%0.0f$'))
                ax[s][l].tick_params(axis='both', which='major', labelsize='larger')
                ax[s][l].tick_params(axis='both', which='minor', labelsize='larger')

                if s == 3:
                    ax[s][l].set_xlabel(r'$\theta$ [arcmin]', size='larger')
                if l == 0:
                    ax[s][l].set_ylabel('%s\n' % self.plotting['redshift_s'][s] + r'Boost factors', size='larger',
                                        linespacing=3)
                if s == 0:
                    ax[s][l].set_title(self.plotting['redshift_l'][l], size='larger')

                ax[s][l].axvspan(self.config['thlims'][0], self.plotting['th_limit'][l], color='gray', alpha=0.2)

        ax[0][len(self.zbins['sbins'])].legend(frameon=False, fontsize=16, loc='lower right')
        self.save_plot('boost_factors')


    def plot_boost_factors_cov(self):
        
        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')
        cmap = self.plotting['cmap']

        fig, ax = plt.subplots(len(self.zbins['sbins']),len(self.zbins['lbins']),figsize = (20, 15), sharey=True,sharex=True) #gridspec_kw = {'height_ratios':[1, 1, 1, 1, 1]})
        fig.subplots_adjust(hspace=0.1,wspace=0.1, right=0.885)

        def corrmatrix(COV):
           COV = np.mat(COV)
           D = np.diag(np.sqrt(np.diag(COV)))
           d = np.linalg.inv(D)
           CORR = d*COV*d
           return CORR

        # To iterate between lens bins
        for l in range(0, len(self.zbins['lbins'])):
           # To iterate between source bins
           for s in range(len(self.zbins['sbins'])):
               # Load cov
               path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
               if os.path.exists(path_test + 'cov_boost_factor'):
                   cov = np.loadtxt(path_test + 'cov_boost_factor')
               
                   corr = corrmatrix(cov)

                   th_edges = np.logspace(np.log10(self.config['thlims'][0]), np.log10(self.config['thlims'][1]), self.config['nthbins']+1)
                   x = th_edges
                   y = th_edges
                   X, Y = np.meshgrid(x, y)
                   ax[s][l].set_xscale('log')
                   ax[s][l].set_yscale('log')
                   im = ax[s][l].pcolormesh(x, y, np.array(corr), vmin = 0.5, vmax=1.)
                   ax[s][l].autoscale('tight')

                   cbar_ax = fig.add_axes([0.9, 0.1, 0.015, 0.8])
                   cbar = fig.colorbar(im, cax=cbar_ax)
                   cbar_ax.set_ylabel(r'$r_{ij}=\mathrm{Corr}(\mathcal{B}_{i}, \mathcal{B}_{j})$', size = 'large')
                   cbar.ax.tick_params(labelsize='large') 
               ax[s][l].set_xlim(self.config['thlims'][0], self.config['thlims'][1])
               ax[s][l].set_ylim(self.config['thlims'][0], self.config['thlims'][1])

               ax[s][l].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%0.0f$'))
               ax[s][l].yaxis.set_major_formatter(ticker.FormatStrFormatter('$%0.0f$'))
               ax[s][l].tick_params(axis='both', which='major', labelsize ='large')
               ax[s][l].tick_params(axis='both', which='minor', labelsize='large')

               if s == 3:
                  ax[s][l].set_xlabel(r'$\theta$ [arcmin]',size='large')
               if l == 0:
                  ax[s][l].set_ylabel('%s\n'%self.plotting['redshift_s'][s] + r'$\theta$ [arcmin]',size='large', linespacing = 3)
               if s == 0:
                  ax[s][l].set_title(self.plotting['redshift_l'][l],size='large')

        self.save_plot('boost_factors_cov')

    def plot_randoms(self):
        """
        Makes plot of the tangential shear around random points.
        """

        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        if self.plotting['use_cmap']:
            cmap = self.plotting['cmap']
            cmap_step = 0.25
            colors = []
            for s in range(len(self.zbins['sbins'])):
                colors.append(plt.get_cmap(cmap)(cmap_step * s))
        else:
            colors = self.plotting['colors']

        
        labels = [self.plotting['catname']]
        c = 0  # If adding im3shape, s=1
        markers = ['o', '^']
        fs = 18  # fontsize
        cmap = self.plotting['cmap']
        cmap_step = 0.25
        fig, ax = plt.subplots(len(self.zbins['sbins']), len(self.zbins['lbins']), figsize=(17.25, 13.8), sharey=True, sharex=True)
        fig.subplots_adjust(hspace=0.0, wspace=0.0)

        for l in range(0, len(self.zbins['lbins'])):
            for s in range(len(self.zbins['sbins'])):
                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                th, gt, err = np.loadtxt(path_test + 'randoms', unpack=True)

                ax[s][l].axhline(y=0, ls=':', color='k')
                ax[s][l].errorbar(th * (1 + 0.07 * c), gt, err, fmt=markers[c], color=colors[s],
                                  mec=colors[s], label=labels[c], capsize=1.3)

                ax[s][l].set_xlim(2.5, 300)
                ax[s][l].set_ylim(-2.3 * 10 ** (-4), 2.3 * 10 ** (-4))
                ax[s][l].set_xscale('log')

                ax[s][l].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%0.0f$'))
                ax[s][l].tick_params(axis='both', which='major', labelsize=fs)
                ax[s][l].tick_params(axis='both', which='minor', labelsize=fs)

                if s == 3:
                    ax[s][l].set_xlabel(r'$\theta$ [arcmin]', fontsize=fs)
                if l == 0:
                    ax[s][l].set_ylabel('%s\n' % self.plotting['redshift_s'][s] + r'$\gamma_t (\theta)$', fontsize=fs,
                                        linespacing=3)
                    ax[s][l].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                if s > 0:
                    ax[s][l].yaxis.get_offset_text().set_visible(False)
                if s == 0:
                    ax[s][l].set_title(self.plotting['redshift_l'][l], fontsize=fs)

        ax[0][len(self.zbins['sbins'])].legend(frameon=False, fancybox=True, prop={'size': 13}, numpoints=1, loc='upper right')
        self.save_plot('plot_randoms')

    def plot_gammax(self):
        """
        Makes plot of the cross-component gammax.
        Top panel: Plot of the gammax measurement for a single lens-source combination.
        Bottom panel: chi2 distribution from all lens-source combinations.
        """

        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        labels = [self.plotting['catname']]
        c = 0  # for metacal
        markers = ['o', '^']
        fs = 12  # fontsize
        cmap = self.plotting['cmap']
        cmap_step = 0.25
        c1 = plt.get_cmap(cmap)(cmap_step * 1.5)
        c2 = plt.get_cmap(cmap)(cmap_step * 3)
        colors = [c1, c2]
        fig, ax = plt.subplots(2, 1, figsize=(5, 8), sharey=False, sharex=False)
        fig.subplots_adjust(hspace=0.2, wspace=0.)

        # TOP panel
        l = 0
        s = 0
        path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
        th, gx, err = np.loadtxt(path_test + 'gx', unpack=True)
        ax[0].axhline(y=0, ls=':', color='k')
        ax[0].errorbar(th * (1 + 0.07 * c), gx * th, err * th, fmt=markers[c], color=colors[c],
                       mec=colors[c], label=labels[c], capsize=1.3)
        ax[0].set_xlim(2.5, 300)
        ax[0].set_ylim(-8 * 10 ** (-3), 8 * 10 ** (-3))
        ax[0].set_xscale('log')

        ax[0].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%0.0f$'))
        ax[0].set_xlabel(r'$\theta$ [arcmin]', fontsize=fs)
        ax[0].set_ylabel(r'$\gamma_{\times} \times \theta$', fontsize=fs)
        ax[0].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax[0].text(0.5, 0.84, self.plotting['redshift_l'][l], horizontalalignment='center',
                   verticalalignment='center', transform=ax[0].transAxes, fontsize=12)
        ax[0].text(0.5, 0.92, self.plotting['redshift_s'][s], horizontalalignment='center',
                   verticalalignment='center', transform=ax[0].transAxes, fontsize=12)
        ax[0].legend(frameon=False, fancybox=True, prop={'size': 10}, numpoints=1, loc='lower left')

        # BOTTOM panel
        chi2s = []
        for l in range(0, len(self.zbins['lbins'])):
            for s in range(len(self.zbins['sbins'])):
                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                chi2, ndf = self.get_chi2(path_test, 'gx')
                chi2s.append(chi2)
        chi2s = np.array(chi2s)
        ax[1].hist(chi2s, bins=7, color=colors[c], ec=colors[c], lw=2, normed=True, histtype='step', alpha=1,
                   label=labels[c])

        # Chi2 distribution
        ndf = len(gx)
        t = np.linspace(0., 60, 300)
        ax[1].fill_between(t, 0, scipy.stats.chi2.pdf(t, ndf), color='gray', alpha=0.5,
                           label='$\chi^2$ pdf (ndf$ = %d$)' % ndf)
        ax[1].set_ylim(0, 0.1)
        ax[1].legend(frameon=False, fancybox=True, prop={'size': 10}, numpoints=1, loc='best')
        ax[1].set_xlabel(r'$\chi^2_\mathrm{Null}$')
        self.save_plot('plot_gammax')


class ResponsesScale(GGL):
    """
    Subclass that obtains the scale dependence responses (NK correlations) for all the lens-source bin combinations. Both for randoms and lenses.
    """

    def __init__(self, basic, config, paths, zbins, plotting):
        GGL.__init__(self, basic, config, paths)
        self.zbins = zbins
        self.plotting = plotting

    def get_path_test_allzbins(self):
        return os.path.join(self.paths['runs_config'], 'responses_nk') + '/'

    def get_path_test(self, lbin, sbin):
        return os.path.join(self.paths['runs_config'], 'responses_nk', lbin + '_' + sbin) + '/'

    def get_path_fidmeasurement(self, lbin, sbin):
        return os.path.join(self.paths['runs_config'], 'measurement', lbin + '_' + sbin) + '/'
    
    def get_twopointfile_name(self, string):
        return os.path.join(self.get_path_test_allzbins() + '%s_nkresp_twopointfile.fits' % string)

    def save_responses_nk(self, path_test, responses, end):
        np.savetxt(path_test + 'responses_nk_no_jackknife_%s' % end, responses,
                   header='theta(arcmin), R_nk, Rgamma_nk, Rs_nk')

    def load_responses_nk(self, path_test, end):
        theta, R_nk, Rgamma_nk, Rs_nk = np.loadtxt(path_test + 'responses_nk_%s' % end, unpack=True)
        return theta, R_nk, Rgamma_nk, Rs_nk

    def load_responses_nk_errors(self, path_test, end):
        theta, R_nk, err = np.loadtxt(path_test + 'mean_R_nk_JK_%s' % end, unpack=True)
        return theta, R_nk, err

    def save_responses_mean(self, responses_mean, end):
        responses_mean = np.array([responses_mean[sbin] for sbin in self.zbins['sbins']])
        np.savetxt(self.get_path_test_allzbins() + 'responses_mean_%s' % end, responses_mean,
                   header='R_mean, Rgamma_mean, Rs_mean')

    def load_responses_mean(self, end):
        R_mean, Rgamma, Rs = np.loadtxt(self.get_path_test_allzbins() + 'responses_mean_%s' % end, unpack=True)
        return R_mean


    def run(self):

        make_directory(self.get_path_test_allzbins())
        for sbin in self.zbins['sbins']:

            if self.basic['mode'] == 'data':
                source_all, source_all_5sels, calibrator = self.load_data_or_sims()
            else:
                raise Exception('It makes no sense to run this test on a simulation. Mode should be data, edit on info.py file.')

            print('Running measurement for source %s.' % sbin)
            print('Source bin:', self.zbins[sbin][0], self.zbins[sbin][1])
            
            source = self.load_metacal_bin(source_all, source_all_5sels, calibrator, bin_low=self.zbins[sbin][0], bin_high=self.zbins[sbin][1], reduce_mem=True, ind_responses=True)
            # R_mean = source['R_mean'] # average for the tomographic bin
            # R = source['R'] # total R for each individual galaxy
            print('Length source', sbin, len(source['ra']))

            # Deleting arrays we do not need before running Treecorr due to memory issues.
            del source_all
            del source_all_5sels
            del calibrator

            # New!! we need to subtract the mean shear # not really needed here but leave it just in case
            source, mean_shear = self.subtract_mean_shear(source)
            for l, lbin in enumerate(self.zbins['lbins']):
                path_test = self.get_path_test(lbin, sbin)
                path_fid = self.get_path_fidmeasurement(lbin, sbin) # where the fiducial measurements are

                if os.path.exists(path_test + 'gt_boosted'):
                    print('Measurements for this bin already exist. SKIPPING!')

                else:
                    print('Running measurement for lens %s.' % lbin)
                    make_directory(path_test)

                    # Preparing arrays for lenses and randoms
                    # ---------------------------------------------
                    lens_all = pf.getdata(self.paths['lens'])
                    random_all = pf.getdata(self.paths['randoms'])
                    lens = lens_all[(lens_all['z'] > self.zbins[lbin][0]) & (lens_all['z'] < self.zbins[lbin][1])]
                    random = random_all[(random_all['z'] > self.zbins[lbin][0]) & (random_all['z'] < self.zbins[lbin][1])]
                    del lens_all
                    del random_all
                    print('Length lens', lbin, len(lens['ra']), self.zbins[lbin][0], self.zbins[lbin][1])
                    print('Length randoms', lbin, len(random['ra']), self.zbins[lbin][0], self.zbins[lbin][1])

                    # Lenses run
                    # ----------------
                    theta, R_nk, weights, npairs = self.run_treecorr_jackknife(lens, source, 'NK_R')
                    theta_fid, gts, gxs, errs, weights_fid, npairs_fid = self.load_runs(path_fid, random_bool=False) # Load from fiducial run
                    assert np.allclose(theta,theta_fid)
                    assert np.allclose(weights, weights_fid)

                    R_jk, _, wnum_jk = self.numerators_jackknife(R_nk, R_nk, weights)
                    Rnum_ex, _, wnum_ex = self.numerators_exact(R_nk, R_nk, weights)

                    # Do this too for the loaded fiducial run to be able to combine here
                    gtnum_jk, gxnum_jk, wnum_jk = self.numerators_jackknife(gts, gxs, weights)
                    gtnum_ex, gxnum_ex, wnum_ex = self.numerators_exact(gts, gxs, weights)

                    # Randoms run
                    # ----------------
                    theta, R_nk_r, weights, npairs = self.run_treecorr_jackknife(random, source, 'NK_R')
                    theta_fid, gts, gxs, errs, weights_fid, npairs_fid = self.load_runs(path_fid, random_bool=True) # Load from fiducial run
                    assert np.allclose(theta,theta_fid)
                    assert np.allclose(weights, weights_fid)

                    R_jk_r, _, wnum_jk_r = self.numerators_jackknife(R_nk_r, R_nk_r, weights)
                    R_ex_r, _, wnum_ex_r = self.numerators_exact(R_nk_r, R_nk_r, weights)

                    gtnum_jk_r, gxnum_jk_r, wnum_jk_r = self.numerators_jackknife(gts, gxs, weights)
                    gtnum_ex_r, gxnum_ex_r, wnum_ex_r = self.numerators_exact(gts, gxs, weights)

                    # Combine final estimator for JK covariance and mean, with and without boost factors
                    # ----------------
                    bf_all = self.compute_boost_factor_jackknife(lens['jk'], random['jk'], wnum_jk, wnum_jk_r, lens['w'])
                    gt_all_boosted = bf_all*gtnum_jk/R_jk - gtnum_jk_r/R_jk_r
                    gx_all = gxnum_jk/R_jk - gxnum_jk_r/R_jk_r

                    # Combine final estimator for exact results (as if done without JK)
                    # ----------------
                    bf = self.compute_boost_factor_exact(wnum_ex, wnum_ex_r, lens, random)
                    gt_boosted = bf*gtnum_ex/Rnum_ex - gtnum_ex_r/R_ex_r 
                    gx = gxnum_ex/Rnum_ex - gxnum_ex_r/R_ex_r

                    # Get the mean and standard deviation of all these quantities and save them in files
                    # ---------------
                    self.process_run(gtnum_ex_r/R_ex_r, gtnum_jk_r/R_jk_r, theta, path_test, 'randoms')
                    self.process_run(bf, bf_all, theta, path_test, 'boost_factor')
                    self.process_run(gt_boosted, gt_all_boosted, theta, path_test, 'gt_boosted')
                    self.process_run(gx, gx_all, theta, path_test, 'gx')
                    self.process_run(Rnum_ex/wnum_ex, R_jk/wnum_jk, theta, path_test, 'R_nk')

    
    def compute_Rs(self, e_ix, delta_gamma):
        """
        Computes R_s.
        e_ix: Dictionary per each component 1p, 1m, 2p, 2m.
        delta_gamma: value of the artificially applied shear to the images. 
        It can be averaged over all angular scales, or averaged in angular bins using NK correlation.
        """
        Rs11_mean = (e_ix['1p'] - e_ix['1m']) / delta_gamma
        Rs22_mean = (e_ix['2p'] - e_ix['2m']) / delta_gamma
        Rs_mean = 0.5 * (Rs11_mean + Rs22_mean)
        return Rs_mean

    def compute_chi2_consistency(self, datavector, constant, cov):
        diff = datavector - constant
        chi2 = np.dot(diff.T, np.dot(np.linalg.inv(cov), diff))

        # Hartlap factor
        N = self.config['njk']  # number of jackknife regions
        p = len(diff)  # number of angular bins
        factor = (N - p - 2) / float(N - 1)
        chi2_hartlap = chi2 * factor
        return chi2_hartlap

    def fit_constant_minuit(self, datavector, cov):
        """
        Function to fit a constant to some data.
        Minimizes chi2 using Minuit
        """
        from iminuit import Minuit
        data_mat = np.mat(datavector)

        def f(c):
            return np.abs(np.array((data_mat - c) * np.linalg.inv(cov) * (data_mat.T - c))[0][0])

        m = Minuit(f, print_level=0, errordef=1, pedantic=False)
        m.migrad()
        fit, err_fit = m.values['c'], m.errors['c']
        hartlap_factor = (self.config['njk'] - len(datavector) - 2) / float(self.config['njk'] - 1)
        chi2_fit = f(fit) * hartlap_factor

        print('fit, err_fit = %0.3e +- %0.3e (Minuit)' % (fit, err_fit))
        print('chi2_fit/ndf: %0.2f/%d (Minuit)' % (chi2_fit, (len(datavector) - 1)))
        return fit, err_fit, chi2_fit, len(datavector) - 1

    def fit_constant_least_squares(self, datavector, cov, R0):
        """
        Function to fit a constant to some data.
        Minimizes chi2 using least squares from scipy
        R0: initial guess for the constant
        """
        from scipy import optimize
        data_mat = np.mat(datavector)

        def f(c):
            return np.abs(np.array((data_mat - c) * np.linalg.inv(cov) * (data_mat.T - c))[0][0])

        # ipdb.set_trace()
        opt = optimize.least_squares(f, R0)
        fit = opt.x
        print('Fit with least squares:', fit)

        hartlap_factor = (self.config['njk'] - len(datavector) - 2) / float(self.config['njk'] - 1)
        chi2_fit = f(fit[0]) * hartlap_factor

        print('fit = %0.3e (least squares)' % (fit[0]))
        print('chi2_fit/ndf: %0.2f/%d (least squares)' % (chi2_fit, (len(datavector) - 1)))
        return fit

    def build_dictionary_e_ix(self, lens, source_sels, average_type):
        """
        Function to build the dictionary e_ix, which is the mean ellipticity for 
        the component i=(1,2) for a given selection x=(p,m). Therefore, it will have
        the mean of e_1p['1p'], e_1m['1m'], e_2p['2p'], e_2m['2m'], which will be obtained
        from the selections 1p, 1m, 2p, 2m respectively. This dictionary is used 
        to compute the selection response Rs.
        """

        e_ix = {}
        components = ['1p', '1m', '2p', '2m']
        for i, comp in enumerate(components):
            source_component_ix = {
                'ra': source_sels['ra'][i],
                'dec': source_sels['dec'][i],
                'e_ix': source_sels['e%s' % comp[0]][
                    i]}  # Choose e1 for 1p, 1m selections, and e2 for 2p, 2m selections.
            if average_type == 'mean':
                e_ix[comp] = np.mean(source_component_ix['e_ix'])
            if average_type == 'NK_no_jackknife':
                theta, e_ix[comp] = self.run_nk_no_jackknife(lens, source_component_ix,
                                                             scalar=source_component_ix['e_ix'])

            if average_type == 'NK_jackknife':
                theta, xi_nk, weights, npairs = self.run_treecorr_jackknife(lens, source_component_ix,
                                                                            type_corr='NK_e_ix')
                e_ixnum, _, wnum = self.numerators_jackknife(xi_nk, xi_nk, weights)
                e_ix[comp] = e_ixnum / wnum  # contains all the jackknife regions measurements

        return e_ix

    def run_responses_tomo(self, lens, source, source_sels, delta_gamma, average_type, path_test, lens_or_random):
        """
        Function that computes mean response for each source bin or scale dependent responses for each lens-source combination.
        Uses NK TreeCorr correlation for scale dependence.
        Source: Source catalog for a certain redshift bin for the unsheared selection. 
        Source_sels: List of source catalogs for each of the four sheared selections (sheared 1p, 1m, 2p, 2m)
                     but with the unsheared quantities ra, dec, e1, e2, to obtain the new selection response.  
        delta_gamma: value of the artificially applied shear to the images. 
        average_type: way to average the ellipticities, can be:
                     - 'mean': np.mean
                     - 'NK_no_jackknife': using NK correlations without jackknife
                     - 'NK_jackknife': using NK correlations with jackknife
        lens_or_random: string for saving.
        Rgamma: (R11 + R22)/2 for each galaxy.
        """

        if average_type == 'mean':
            Rgamma = np.mean(source['Rgamma'])

        if average_type == 'NK_no_jackknife':
            theta, Rgamma = self.run_nk_no_jackknife(lens, source, scalar=source['Rgamma'])

        if average_type == 'NK_jackknife':
            theta, xi_nk, weights, npairs = self.run_treecorr_jackknife(lens, source, type_corr='NK_Rgamma')
            Rgammanum, _, wnum = self.numerators_jackknife(xi_nk, xi_nk, weights)
            Rgamma = Rgammanum / wnum  # contains all the jackknife regions measurements (i.e. Rgamma_all)

        e_ix = self.build_dictionary_e_ix(lens, source_sels, average_type)
        Rs = self.compute_Rs(e_ix, delta_gamma)

        R = Rgamma + Rs

        print('average_type, R, Rgamma, Rs', average_type, R, Rgamma, Rs)

        if average_type == 'mean':
            responses = [R, Rgamma, Rs]
            return responses

        if average_type == 'NK_no_jackknife':
            responses = zip(theta, R, Rgamma, Rs)
            self.save_responses_nk(path_test, responses, lens_or_random)

        if average_type == 'NK_jackknife':
            self.process_run(R, theta, path_test, 'R_nk_JK_%s' % lens_or_random)
            self.process_run(Rgamma, theta, path_test, 'Rgamma_nk_JK_%s' % lens_or_random)
            self.process_run(Rs, theta, path_test, 'Rs_nk_JK_%s' % lens_or_random)

    def run_deprecated(self):
        """
        Runs the NK responses between lenses and sources.
        Runs for lenses and randoms too.
        Uses old (explicit) functions to compute the selection response.
        """
        lens_all, random_all, source_all, source_all_5sels, calibrator = self.load_data_or_sims()
        resp = {}
        responses_mean = {}
        for sbin in self.zbins['sbins']:

            print('Running responses test for source %s.' % sbin)
            source = self.load_metacal_bin(source_all, source_all_5sels, calibrator, bin_low=self.zbins[sbin][0],
                                           bin_high=self.zbins[sbin][1])
            resp[sbin] = source['Rmean']
            print('R = source[Rmean]', resp)
            source_sels = self.load_metacal_bin_sels_responses(source_all_5sels, bin_low=self.zbins[sbin][0],
                                                               bin_high=self.zbins[sbin][1])
            delta_gamma = 2 * 0.01

            for lbin in self.zbins['lbins']:
                print('Running responses test for lens %s.' % lbin)
                path_test = self.get_path_test(lbin, sbin)
                make_directory(path_test)

                lens = lens_all[(lens_all['z'] > self.zbins[lbin][0]) & (lens_all['z'] < self.zbins[lbin][1])]
                responses_mean[sbin] = self.run_responses_tomo(lens, source, source_sels, delta_gamma,
                                                               average_type='mean', path_test=path_test,
                                                               lens_or_random='lens')
                # comment:responses_nk_no_jackknife = self.run_responses_tomo(lens, source, source_sels, delta_gamma, average_type='NK_no_jackknife', path_test=path_test, lens_or_random='lens') #much slower
                self.run_responses_tomo(lens, source, source_sels, delta_gamma, average_type='NK_jackknife',
                                        path_test=path_test, lens_or_random='lens')

                # comment:random = random_all[(random_all['z'] > self.zbins[lbin][0]) & (random_all['z'] < self.zbins[lbin][1])]
                # comment:responses_nk_no_jackknife = self.run_responses_tomo(random, source, source_sels, delta_gamma, average_type='NK_no_jackknife', path_test=path_test, lens_or_random='random') #super slow
                # comment:self.run_responses_tomo(random, source, source_sels, delta_gamma, average_type='NK_jackknife', path_test=path_test, lens_or_random='random')

        print(resp)
        print(responses_mean)
        self.save_responses_mean(resp, 'destest')
        self.save_responses_mean(responses_mean, 'xcorr')

    def plot(self, lens_random, mask_scales):
        """
        Makes plot comparing the NK responses to the mean ones.
        lens_random: string, can be lens or random.
        mask_scales: boolean, scales to use for the chi2
        Indicates which is the foreground sample when computing the NK correlations.
        """

        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        cmap = self.plotting['cmap']
        cmap_step = 0.25
        c1 = plt.get_cmap(cmap)(0.)
        c2 = plt.get_cmap(cmap)(0.25)
        c3 = plt.get_cmap(cmap)(0.5)
        fig, ax = plt.subplots(len(self.zbins['sbins']), len(self.zbins['lbins']), figsize=(16.5, 13.2), sharey='row',
                               sharex=True)
        fig.subplots_adjust(hspace=0.1, wspace=0.1)

        R_mean_all = self.load_responses_mean('xcorr')

        err_Rmean = np.array([0.0003528, 0.0004943, 0.0004696,
                              0.0005272])  # values from Marco from JK, see slack on 14th of March, propagate errors
        print(R_mean_all)
        for l in range(len(self.zbins['lbins'])):

            for s in range(len(self.zbins['sbins'])):
                R_mean = R_mean_all[s]
                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                # theta, R_nk, _, _ = self.load_responses_nk(path_test, lens_random)
                theta, R_nk_jk, err = self.load_responses_nk_errors(path_test, lens_random)
                cov = np.loadtxt(path_test + 'cov_R_nk_JK_%s' % lens_random)

                if mask_scales:
                    ax[s][l].axvspan(self.config['thlims'][0] * 0.8, self.plotting['th_limit'][l], color='gray',
                                     alpha=0.2)
                    mask = theta > self.plotting['th_limit'][l]
                    chi2 = self.compute_chi2_consistency(datavector=R_nk_jk[mask], constant=R_mean,
                                                         cov=(cov[mask].T)[mask].T)
                    save = 'mask_scales'
                    ndf = len(R_nk_jk[mask])
                    c, err_c, chi2_c, ndf_c = self.fit_constant_minuit(datavector=R_nk_jk[mask],
                                                                       cov=(cov[mask].T)[mask].T)

                else:
                    chi2 = self.compute_chi2_consistency(datavector=R_nk_jk, constant=R_mean, cov=cov)
                    save = 'all_scales'
                    ndf = len(R_nk_jk)
                    c, err_c, chi2_c, ndf_c = self.fit_constant_minuit(datavector=R_nk_jk, cov=cov)
                    c_ls = self.fit_constant_least_squares(datavector=R_nk_jk, cov=cov, R0=R_mean)

                ax[s][l].margins(x=0, y=0.3)
                t = np.linspace(0, max(theta) * 2, 10)
                ax[s][l].plot(t, [R_mean] * len(t), '-', lw=2, color=c1, mec=c1, label=r'$R_{\mathrm{mean}}$')
                ax[s][l].fill_between(t, np.array([R_mean - err_Rmean[s] for i in t]),
                                      np.array([R_mean + err_Rmean[s] for i in t]), alpha=0.4, edgecolor=c1,
                                      facecolor=c1)
                # ax[s][l].plot(theta, R_nk, '-', lw=2, color=c2, mec=c2, label=r'$R_{\mathrm{nk}}$')
                # ax[s][l].plot(theta, [np.mean(R_nk)] * len(theta), '--', lw=2, color=c2, mec=c2,
                #              label=r'$\overline{R_{\mathrm{nk}}}$')
                ax[s][l].errorbar(theta, R_nk_jk, err, fmt='o', color=c3, mec=c3, markersize=3.,
                                  label=r'$R_{\mathrm{nk, jk}}$')
                # ax[s][l].plot(theta, [np.mean(R_nk_jk)] * len(theta), ':', lw=2, color=c3, mec=c3,
                #              label=r'$\overline{R_{\mathrm{nk, jk}}}$')

                # ax[s][l].plot(theta, [c_ls] * len(theta), '--', lw=2, color=c3, mec=c3,
                #              label=r'Fit with least squares')
                ax[s][l].fill_between(t, np.array([c - err_c for i in t]),
                                      np.array([c + err_c for i in t]), alpha=0.4, edgecolor=c3, facecolor=c3,
                                      label='Fit to a constant')

                ax[s][l].set_xscale('log')
                ax[s][l].set_xlim(self.config['thlims'][0], self.config['thlims'][1])
                ax[s][l].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%0.0f$'))
                ax[s][l].tick_params(axis='both', which='major', labelsize='larger')
                ax[s][l].tick_params(axis='both', which='minor', labelsize='larger')

                if s == 3:
                    ax[s][l].set_xlabel(r'$\theta$ [arcmin]', size='larger')
                if l == 0:
                    ax[s][l].set_ylabel('%s\n' % self.plotting['redshift_s'][s] + r'Responses', size='larger',
                                        linespacing=3)
                if s == 0:
                    ax[s][l].set_title(self.plotting['redshift_l'][l], size='larger')

                diff = R_mean / R_nk_jk - 1

                # ax[s][l].text(0.5, 0.88,
                #              r'Mean $R_{\mathrm{mean}}/R_{\mathrm{nk}}-1 = %0.2f \%%$' % (100 * np.mean(diff)),
                #              horizontalalignment='center', verticalalignment='center', transform=ax[s][l].transAxes,
                #              fontsize='medium')

                # ax[s][l].text(0.5, 0.79,
                #              r'Max $R_{\mathrm{mean}}/R_{\mathrm{nk}}-1 = %0.2f \%%$' % (100 * np.max(np.absolute(diff))),
                #              horizontalalignment='center', verticalalignment='center', transform=ax[s][l].transAxes,
                #              fontsize='medium')

                ax[s][l].text(0.5, 0.85,
                              r'$\chi^2$/ndf$ = %0.1f/%d$' % (chi2, ndf),
                              horizontalalignment='center', verticalalignment='center', transform=ax[s][l].transAxes,
                              fontsize='medium')

        ax[0][4].legend(frameon=False, fontsize=10, loc='lower right')
        self.save_plot('plot_responses_scale_dependence_%s_%s' % (save, lens_random))

    def plot_sigmas(self, lens_random, mask_scales, string):
        """
        Makes plot comparing the NK responses to the mean ones, divided by the uncertainty on the measurement. 
        lens_random: string, can be lens or random.
        string: indicates which twopointfile to load the sigmas from, i.e. from gt, from gt boosted etc.
        Indicates which is the foreground sample when computing the NK correlations.
        """

        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        cmap = self.plotting['cmap']
        cmap_step = 0.25
        c1 = plt.get_cmap(cmap)(0.)
        c2 = plt.get_cmap(cmap)(0.6)
        fig, ax = plt.subplots(len(self.zbins['sbins']), len(self.zbins['lbins']), figsize=(16.5, 13.2), sharey=True,
                               sharex=True)
        fig.subplots_adjust(hspace=0.1, wspace=0.1)

        R_mean_all = self.load_responses_mean('xcorr')
        print(R_mean_all)

        measurement = Measurement(self.basic, self.config, self.paths, self.zbins, self.plotting)
        gammat_file = measurement.load_twopointfile(string)
        gammat = gammat_file.spectra[0]

        for l in range(len(self.zbins['lbins'])):
            for s in range(len(self.zbins['sbins'])):
                R_mean = R_mean_all[s]
                path_test = self.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                theta, R_nk, _ = self.load_responses_nk_errors(path_test, lens_random)
                path_test_measurement = measurement.get_path_test(self.zbins['lbins'][l], self.zbins['sbins'][s])
                cov = np.loadtxt(path_test_measurement + 'cov_gt')
                _, gt = gammat.get_pair(l + 1, s + 1)
                err = gammat.get_error(l + 1, s + 1)
                diff_err = (R_mean / R_nk - 1) * gt / err
                diff = (R_mean / R_nk - 1) * gt

                if mask_scales:
                    ax[s][l].axvspan(self.config['thlims'][0] * 0.8, self.plotting['th_limit'][l], color='gray',
                                     alpha=0.2)
                    mask = theta > self.plotting['th_limit'][l]
                    chi2 = self.compute_chi2_consistency(datavector=R_mean / R_nk[mask] * gt[mask], constant=gt[mask],
                                                         cov=(cov[mask].T)[mask].T)
                    save = 'mask_scales'
                    ndf = len(R_nk[mask])
                else:
                    chi2 = self.compute_chi2_consistency(datavector=R_mean / R_nk * gt, constant=gt, cov=cov)
                    save = 'all_scales'
                    ndf = len(R_nk)

                print('l, s, chi2:', l, s, chi2)
                ax[s][l].plot(theta, diff_err, lw=2, color=c1, mec=c1)
                ax[s][l].axhline(y=0, color='k', ls=':')
                # ax[s][l].plot(theta, err, lw=2, color=c2, mec=c2, label=r'$\sigma_{\gamma_t, \mathrm{JK}}$')
                # ax[s][l].plot(theta, diff/err, lw=2, color=c2, mec=c2, label=r'$(R_{\mathrm{nk}} - R_{\mathrm{mean}})/\sigma_{\gamma_t}$')
                ax[s][l].set_xscale('log')
                # ax[s][l].set_yscale('log')
                ax[s][l].set_xlim(self.config['thlims'][0], self.config['thlims'][1])
                ax[s][l].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%0.0f$'))
                ax[s][l].tick_params(axis='both', which='major', labelsize='larger')
                ax[s][l].tick_params(axis='both', which='minor', labelsize='larger')

                if s == 3:
                    ax[s][l].set_xlabel(r'$\theta$ [arcmin]', size='larger')
                if l == 0:
                    ax[s][l].set_ylabel('%s\n' % self.plotting['redshift_s'][
                        s] + r'$(R_{\mathrm{mean}}/R_{\mathrm{nk}}-1)\gamma_t/\sigma_{\gamma_t, \mathrm{JK}}$',
                                        size='larger', linespacing=3)
                if s == 0:
                    ax[s][l].set_title(self.plotting['redshift_l'][l], size='larger')

                ax[s][l].text(0.5, 0.87,
                              r'$\chi^2$/ndf$ = %0.3f/%d$' % (chi2, ndf),
                              horizontalalignment='center', verticalalignment='center', transform=ax[s][l].transAxes,
                              fontsize='medium')

        ax[0][4].legend(frameon=False, fontsize=16, loc='lower right')
        self.save_plot('plot_responses_%s_diff_%s' % (lens_random, string))


class TestStars(GGL):
    """
    SubClass to test if the tangential shear around stars is consistent with zero.
    Uses no tomography for the source sample.
    """

    def __init__(self, basic, config, paths, zbins, plotting):
        GGL.__init__(self, basic, config, paths)
        self.zbins = zbins
        self.plotting = plotting

    def get_path_test(self, typestars):
        return os.path.join(self.paths['runs_config'], 'stars_%s' % typestars) + '/'

    def run(self, typestars):

        """
        Runs the gglensing measurement with stars as lenses and the full source sample.
        typestars: string that indicates which stars to load. Either 'bright' or 'faint'.
        Includes:
        - Mean response calculation for the full sample.
        - Random points subtraction.
        - Jackknife covariance calculation.
        """

        stars = pf.getdata(self.paths['y1'] + 'star_%s.fits' % typestars)
        randoms = pf.getdata(self.paths['y1'] + 'random.fits')
        mask_randoms = np.random.randint(len(randoms), size=len(stars) * 10)
        sources = pf.getdata(self.paths['y1'] + 'metacal_sel_allbins.fits')

        path_test = self.get_path_test(typestars)
        make_directory(path_test)

        theta, gts, gxs, errs, weights, npairs = self.run_treecorr_jackknife(stars, sources, 'NG')
        gtnum, gxnum, wnum = self.numerators_jackknife(gts, gxs, weights)

        theta, gts, gxs, errs, weights, npairs = self.run_treecorr_jackknife(randoms[mask_randoms], sources, 'NG')
        gtnum_r, gxnum_r, wnum_r = self.numerators_jackknife(gts, gxs, weights)

        R = self.run_responses_mean_notomo(sources['Rgamma'])

        gt_all = (gtnum / wnum) / R - (gtnum_r / wnum_r) / R
        gx_all = (gxnum / wnum) / R - (gxnum_r / wnum_r) / R

        self.process_run(gt_all, theta, path_test, 'gt')
        self.process_run(gx_all, theta, path_test, 'gx')

    def plot(self):
        """
        Make plot for all the stars test.
        """

        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        types = ['bright', 'faint']
        styles = ['o', '^']
        shifts = [0, 1, 2, 3]
        colors = [plt.get_cmap(self.plotting['cmap'])(0), plt.get_cmap(self.plotting['cmap'])(0.7)]
        titles_c = [self.plotting['catname']]
        s = 0  # If adding im3shape, s=1

        fig, ax = plt.subplots(1, 1, figsize=(5., 5.), sharey=False, sharex=False)
        ax.axhline(y=0, ls=':', color='k')

        for t, type in enumerate(types):

            path_test = self.get_path_test(type)
            th, gt, err = np.loadtxt(path_test + 'mean_gt', unpack=True)

            if s == 0:
                ax.errorbar(th * (1 + 0.05 * (shifts[t])), gt, err, fmt=styles[t], color=colors[s], mec=colors[s],
                            label=titles_c[s] + ' ' + type, markersize=4.5)
            if s == 1:
                ax.errorbar(th * (1 + 0.05 * (shifts[t + 2])), gt, err, fmt=styles[t], color=colors[s], mec=colors[s],
                            label=titles_c[s] + ' ' + type, markersize=4.5)

            chi2, ndf = self.get_chi2(path_test, 'gt')

        ax.set_xlim(2.5, 300)
        ax.set_ylim(-2 * 10 ** (-4), 2 * 10 ** (-4))
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%d$'))
        ax.set_xlabel(r'$\theta$ [arcmin]', size='large')
        ax.set_ylabel(r'$\gamma_{t} (\theta)$', size='large')
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax.legend(frameon=False, fancybox=True, prop={'size': 12}, numpoints=1, loc='best')
        self.save_plot('plot_stars')


class TestPSF():
    """
    SubClass to test if the psf residuals are compatible with zero.
    Uses no tomography for the lens sample.
    """

    def __init__(self, ggl_class, config, paths, zbins, plotting):
        self.zbins = zbins
        self.plotting = plotting
        self.ggl = ggl_class
        self.paths = paths
        self.config = config

    def get_path_test(self, band):
        return os.path.join(self.paths['runs_config'], 'psfresiduals_band_', band) + '/'

    def get_path_test_Tres(self, band):
        return os.path.join(self.paths['runs_config'], 'Tresiduals_band_', band) + '/'

    def save_psf_residuals_y3(self, ra_lims, dec_lims):
        """
        Computes the psf residuals for the r band and saves them to a file.
        """

        min_ra, max_ra = ra_lims
        min_dec, max_dec = dec_lims
        ra_all_bandr, dec_all_bandr, psf1_all_bandr, psf2_all_bandr, res1_all_bandr, res2_all_bandr, mag_all_bandr, exp_all_bandr, band_all_bandr, T_model_bandr, dT_model_bandr = [], [], [], [], [], [], [], [], [], [], []

        ra_all_bandi, dec_all_bandi, psf1_all_bandi, psf2_all_bandi, res1_all_bandi, res2_all_bandi, mag_all_bandi, exp_all_bandi, band_all_bandi, T_model_bandi, dT_model_bandi = [], [], [], [], [], [], [], [], [], [], []

        ra_all_bandz, dec_all_bandz, psf1_all_bandz, psf2_all_bandz, res1_all_bandz, res2_all_bandz, mag_all_bandz, exp_all_bandz, band_all_bandz, T_model_bandz, dT_model_bandz = [], [], [], [], [], [], [], [], [], [], []

        exposures = []

        y3_exp_dir = self.paths['y3_exp']

        for root, directories, filenames in os.walk(y3_exp_dir):
            for directory in directories:
                exposures.append(directory)

        exposures_int = np.array([int(exp) for exp in exposures])
        exp_unique = np.unique(exposures_int)

        print(len(exp_unique))

        for i in range(len(exp_unique)):

            if np.mod(i, 1000) == 0:
                print(i)

            filename_h = y3_exp_dir + str(exp_unique[i]) + '/exp_psf_cat_%s.fits' % exp_unique[i]
            read_file = fio.FITS(filename_h)
            data = read_file[1]
            info = read_file[2]
            flag = info['flag'].read()

            if np.sum(flag) == 0:
                band = (info['band'].read())[0][0]
                ra = data['ra'].read()
                dec = data['dec'].read()
                obs_e1 = data['obs_e1'].read()
                obs_e2 = data['obs_e2'].read()
                piff_e1 = data['piff_e1'].read()
                piff_e2 = data['piff_e2'].read()
                obs_flag = data['obs_flag'].read()

                piff_T = data['piff_T'].read()
                obs_T = data['obs_T'].read()

                res_T = piff_T - obs_T

                mag = data['mag'].read()

                # Resiual psf is the difference between the measurement of the psf (e1) and the model of the psf(psfex)
                # at the position of the stars, that is the only place you can measure the psf
                res1 = obs_e1 - piff_e1
                res2 = obs_e2 - piff_e2

                # Use only reserved stars and spt
                mask = ((dec > min_dec) & (dec < max_dec) & (obs_flag == 64))

                if band == 'r':
                    ra_all_bandr.extend(ra[mask])
                    dec_all_bandr.extend(dec[mask])
                    psf1_all_bandr.extend(piff_e1[mask])
                    psf2_all_bandr.extend(piff_e2[mask])
                    res1_all_bandr.extend(res1[mask])
                    res2_all_bandr.extend(res2[mask])
                    mag_all_bandr.extend(mag[mask])
                    T_model_bandr.extend(piff_T[mask])
                    dT_model_bandr.extend(res_T[mask])

                elif band == 'i':
                    ra_all_bandi.extend(ra[mask])
                    dec_all_bandi.extend(dec[mask])
                    psf1_all_bandi.extend(piff_e1[mask])
                    psf2_all_bandi.extend(piff_e2[mask])
                    res1_all_bandi.extend(res1[mask])
                    res2_all_bandi.extend(res2[mask])
                    mag_all_bandi.extend(mag[mask])
                    T_model_bandi.extend(piff_T[mask])
                    dT_model_bandi.extend(res_T[mask])


                elif band == 'z':
                    ra_all_bandz.extend(ra[mask])
                    dec_all_bandz.extend(dec[mask])
                    psf1_all_bandz.extend(piff_e1[mask])
                    psf2_all_bandz.extend(piff_e2[mask])
                    res1_all_bandz.extend(res1[mask])
                    res2_all_bandz.extend(res2[mask])
                    mag_all_bandz.extend(mag[mask])
                    T_model_bandz.extend(piff_T[mask])
                    dT_model_bandz.extend(res_T[mask])

                else:
                    print('no correct band alloted')

        # Convert to numpy array and simplify name
        # plt.hist(mag_all, bins=25)
        ra_bandr = np.array(ra_all_bandr)
        dec_bandr = np.array(dec_all_bandr)
        psf1_bandr = np.array(psf1_all_bandr)
        psf2_bandr = np.array(psf2_all_bandr)
        res1_bandr = np.array(res1_all_bandr)
        res2_bandr = np.array(res2_all_bandr)

        print('Number of exposures:', len(ra_bandr))

        c1 = pf.Column(name='RA', format='E', array=ra_bandr)
        c2 = pf.Column(name='DEC', format='E', array=dec_bandr)
        c3 = pf.Column(name='E1', format='E', array=res1_bandr)
        c4 = pf.Column(name='E2', format='E', array=res2_bandr)

        CC = [c1, c2, c3, c4]
        hdu = pf.BinTableHDU.from_columns(CC, nrows=len(ra_bandr))
        hdu.writeto(self.paths['y3'] + 'psfresiduals_bandr.fits', clobber=True)

        c1 = pf.Column(name='RA', format='E', array=ra_bandr)
        c2 = pf.Column(name='DEC', format='E', array=dec_bandr)
        c3 = pf.Column(name='E1', format='E',
                       array=np.array(psf1_all_bandr) * np.array(dT_model_bandr) / np.array(T_model_bandr))
        c4 = pf.Column(name='E2', format='E',
                       array=np.array(psf2_all_bandr) * np.array(dT_model_bandr) / np.array(T_model_bandr))

        CC = [c1, c2, c3, c4]
        hdu = pf.BinTableHDU.from_columns(CC, nrows=len(ra_bandr))
        hdu.writeto(self.paths['y3'] + 'Tresiduals_bandr.fits', clobber=True)

        ra_bandi = np.array(ra_all_bandi)
        dec_bandi = np.array(dec_all_bandi)
        psf1_bandi = np.array(psf1_all_bandi)
        psf2_bandi = np.array(psf2_all_bandi)
        res1_bandi = np.array(res1_all_bandi)
        res2_bandi = np.array(res2_all_bandi)

        print('Number of exposures bandi:', len(ra_bandi))

        c1 = pf.Column(name='RA', format='E', array=ra_bandi)
        c2 = pf.Column(name='DEC', format='E', array=dec_bandi)
        c3 = pf.Column(name='E1', format='E', array=res1_bandi)
        c4 = pf.Column(name='E2', format='E', array=res2_bandi)

        CC = [c1, c2, c3, c4]
        hdu = pf.BinTableHDU.from_columns(CC, nrows=len(ra_bandi))
        hdu.writeto(self.paths['y3'] + 'psfresiduals_bandi.fits', clobber=True)

        c1 = pf.Column(name='RA', format='E', array=ra_bandi)
        c2 = pf.Column(name='DEC', format='E', array=dec_bandi)
        c3 = pf.Column(name='E1', format='E',
                       array=np.array(psf1_all_bandi) * np.array(dT_model_bandi) / np.array(T_model_bandi))
        c4 = pf.Column(name='E2', format='E',
                       array=np.array(psf2_all_bandi) * np.array(dT_model_bandi) / np.array(T_model_bandi))

        CC = [c1, c2, c3, c4]
        hdu = pf.BinTableHDU.from_columns(CC, nrows=len(ra_bandi))
        hdu.writeto(self.paths['y3'] + 'Tresiduals_bandi.fits', clobber=True)

        ra_bandz = np.array(ra_all_bandz)
        dec_bandz = np.array(dec_all_bandz)
        psf1_bandz = np.array(psf1_all_bandz)
        psf2_bandz = np.array(psf2_all_bandz)
        res1_bandz = np.array(res1_all_bandz)
        res2_bandz = np.array(res2_all_bandz)

        # w = np.ones(len(ra))
        print('Number of exposures bandz:', len(ra_bandz))

        c1 = pf.Column(name='RA', format='E', array=ra_bandz)
        c2 = pf.Column(name='DEC', format='E', array=dec_bandz)
        c3 = pf.Column(name='E1', format='E', array=res1_bandz)
        c4 = pf.Column(name='E2', format='E', array=res2_bandz)

        CC = [c1, c2, c3, c4]
        hdu = pf.BinTableHDU.from_columns(CC, nrows=len(ra_bandz))
        hdu.writeto(self.paths['y3'] + 'psfresiduals_bandz.fits', clobber=True)

        c1 = pf.Column(name='RA', format='E', array=ra_bandz)
        c2 = pf.Column(name='DEC', format='E', array=dec_bandz)
        c3 = pf.Column(name='E1', format='E',
                       array=np.array(psf1_all_bandz) * np.array(dT_model_bandz) / np.array(T_model_bandz))
        c4 = pf.Column(name='E2', format='E',
                       array=np.array(psf2_all_bandz) * np.array(dT_model_bandz) / np.array(T_model_bandz))

        CC = [c1, c2, c3, c4]
        hdu = pf.BinTableHDU.from_columns(CC, nrows=len(ra_bandz))
        hdu.writeto(self.paths['y3'] + 'Tresiduals_bandz.fits', clobber=True)

    def run_y3(self, bands):
        """
        Obtains the tangential component of the psf residuals around lenses, with random point subtraction.
        Obtains the corresponding jackknife covariance.
        """

        lens = self.ggl.lens
        masklens = ((lens['z'] > zbins['l1'][0]) & (lens['z'] < zbins['l5'][1]))
        lens = lens[masklens]

        random = self.ggl.random
        maskrandom = ((random['z'] > zbins['l1'][0]) & (random['z'] < zbins['l5'][1]))
        print(zbins['l1'][0])
        print(zbins['l5'][1])
        random = random[maskrandom]

        for band in bands:
            print(band)

            psfres = pf.open(self.paths['y3'] + 'psfresiduals_band' + band + '.fits')[1].data

            path_test = self.get_path_test(band)
            make_directory(path_test)

            print('PSF residuals around lenses...')
            # ipdb.set_trace()
            theta, gts, gxs, errs, weights, npairs = self.ggl.run_treecorr_jackknife(lens, psfres, 'NG')
            gtnum, gxnum, wnum = self.ggl.numerators_jackknife(gts, gxs, weights)

            print('PSF residuals around randoms...')
            theta, gts, gxs, errs, weights, npairs = self.ggl.run_treecorr_jackknife(random, psfres, 'NG')
            gtnum_r, gxnum_r, wnum_r = self.ggl.numerators_jackknife(gts, gxs, weights)

            gt_all = gtnum / wnum - gtnum_r / wnum_r

            self.ggl.process_run(gt_all, theta, path_test, 'gt')

    def run_y3_Tres(self, bands):
        """
        Obtains the tangential component of the psf residuals around lenses, with random point subtraction.
        Obtains the corresponding jackknife covariance.
        """

        lens = self.ggl.lens
        masklens = ((lens['z'] > zbins['l1'][0]) & (lens['z'] < zbins['l5'][1]))
        lens = lens[masklens]

        random = self.ggl.random
        maskrandom = ((random['z'] > zbins['l1'][0]) & (random['z'] < zbins['l5'][1]))
        print(zbins['l1'][0])
        print(zbins['l5'][1])
        random = random[maskrandom]

        for band in bands:
            print(band)

            psfres = pf.open(self.paths['y3'] + 'Tresiduals_band' + band + '.fits')[1].data

            path_test = self.get_path_test_Tres(band)
            make_directory(path_test)

            print('T residuals around lenses...')
            # ipdb.set_trace()
            theta, gts, gxs, errs, weights, npairs = self.ggl.run_treecorr_jackknife(lens, psfres, 'NG')
            gtnum, gxnum, wnum = self.ggl.numerators_jackknife(gts, gxs, weights)

            print('T residuals around randoms...')
            theta, gts, gxs, errs, weights, npairs = self.ggl.run_treecorr_jackknife(random, psfres, 'NG')
            gtnum_r, gxnum_r, wnum_r = self.ggl.numerators_jackknife(gts, gxs, weights)

            gt_all = gtnum / wnum - gtnum_r / wnum_r

            self.ggl.process_run(gt_all, theta, path_test, 'gt')

    def plot(self, bands):
        """
        Makes plot of the psf resdiuals.
        """

        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        fig.subplots_adjust(hspace=0.0, wspace=0.00)

        colors = ['red', 'blue', 'black']
        k = 0

        ax.text(0.7, 0.19, r'Null $\chi^2$/ndf ',
                horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=10)

        for band in bands:
            cmap = plotting['cmap']
            c1 = plt.get_cmap(cmap)(0)
            titles_l = r'$0.15 < z < 0.90 $'
            title_redmagic = 'redMaGiC'

            path_test = self.get_path_test(band)
            th, gt, err = np.loadtxt(path_test + 'mean_gt', unpack=True)
            ax.errorbar((1.025 ** k) * th, gt, err, fmt='.', color=colors[k], mec=c1, markersize=5.7, capsize=1.4,
                        label='Band ' + band)
            chi2, ndf = self.ggl.get_chi2(path_test, 'gt')
            ax.text(0.7, 0.14 - 0.055 * k, r'Band ' + band + ': ' '$%0.1f/%d$' % (chi2, ndf),
                    horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=9,
                    color=c1)
            k += 1
        ax.set_xlim(2.5, 250)
        ax.set_ylim(-1.8 * 10 ** (-5), 1.8 * 10 ** (-5))
        ax.set_xscale('log')
        ax.axhline(y=0, ls=':', color='k')
        ax.text(0.5, 0.85, titles_l, horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(0.5, 0.92, title_redmagic, horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, fontsize=12)

        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%d$'))

        ax.set_xlabel(r'$\theta$ [arcmin]', fontsize=16)
        ax.set_ylabel(r'$\gamma_{t,\mathrm{PSF\, residuals}}$', fontsize=16)

        ax.legend(loc='upper right', fontsize=14, frameon=False)

        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        self.ggl.save_plot('plot_psfresiduals_bandall')

    def plot_Tres(self, bands):
        """
        Makes plot of the psf resdiuals.
        """

        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        fig.subplots_adjust(hspace=0.0, wspace=0.00)

        colors = ['red', 'blue', 'black']
        k = 0

        ax.text(0.7, 0.19, r'Null $\chi^2$/ndf ',
                horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=10)

        for band in bands:
            cmap = plotting['cmap']
            c1 = plt.get_cmap(cmap)(0)
            titles_l = r'$0.15 < z < 0.90 $'
            title_redmagic = 'redMaGiC'

            path_test = self.get_path_test_Tres(band)
            th, gt, err = np.loadtxt(path_test + 'mean_gt', unpack=True)
            ax.errorbar((1.025 ** k) * th, gt, err, fmt='.', color=colors[k], mec=c1, markersize=5.7, capsize=1.4,
                        label='Band ' + band)
            chi2, ndf = self.ggl.get_chi2(path_test, 'gt')
            ax.text(0.7, 0.14 - 0.055 * k, r'Band ' + band + ': ' '$%0.1f/%d$' % (chi2, ndf),
                    horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=9,
                    color=c1)
            k += 1
        ax.set_xlim(2.5, 250)
        ax.set_ylim(-1.8 * 10 ** (-5), 1.8 * 10 ** (-5))
        ax.set_xscale('log')
        ax.axhline(y=0, ls=':', color='k')
        ax.text(0.5, 0.85, titles_l, horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(0.5, 0.92, title_redmagic, horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, fontsize=12)

        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%d$'))

        ax.set_xlabel(r'$\theta$ [arcmin]', fontsize=16)
        ax.set_ylabel(r'$\gamma_{t,\mathrm{T\, residuals}}$', fontsize=16)

        ax.legend(loc='upper right', fontsize=14, frameon=False)

        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        self.ggl.save_plot('plot_Tresiduals_bandall')


class TestSizeSNR(GGL):
    """
    SubClass to test if the tangential shear signal has no dependence on source size or S/N. Using the first lens bin and all sources, as in Y1.
    The variable size_snr on run will be 'size' or 'snr' and will specify the test.
    """

    def __init__(self, basic, config, paths, zbins, plotting, source_nofz_pars):
        GGL.__init__(self, basic, config, paths)
        self.zbins = zbins
        self.plotting = plotting
        self.source_nofz_pars = source_nofz_pars

    def get_path_test(self, size_snr):
        return os.path.join(self.paths['runs_config'], '%s' % size_snr) + '/'

    def save_size_snr(self, path_test, size_snr, result_data, result_theory):
        np.savetxt(path_test + '%s_data' % size_snr, result_data, header='ratio err_ratio')
        np.savetxt(path_test + '%s_theory' % size_snr, result_theory, header='ratio err_ratio')

    def load_size_snr(self, path_test, size_snr):
        ratio_data, err_ratio_data = np.loadtxt(path_test + '%s_data' % size_snr, unpack=True)
        ratio_theory, err_ratio_theory = np.loadtxt(path_test + '%s_theory' % size_snr, unpack=True)
        return ratio_data, err_ratio_data, ratio_theory, err_ratio_theory

    def load_metacal_bin_size_snr(self, source, source_5sels, calibrator, size_snr, high_low, bin_low, bin_high):
        """
        source: dictionary containing relevant columns for the sources, with the baseline selection applied already.
        source_5sels: dictionary containing relevant columns for the sources, with the baseline selection applied already,
                     for each of the 5 selections 1p, 1m, 2p, 2m. 
        calibrator: class to compute the response. Taken from baseline selection.
        zlim_low, zlim_high: limits to select the tomographic bin.
        Obtains 5 masks (unsheared, sheared 1p, 1m, 2p, 2m) to obtain the new selection response.
        Returns: Source dictionary masked with the unsheared mask and with the mean response updated.
        """
        if high_low == 'low':
                masks = [(source_5sels['sheared']['som_bin'][i] >= bin_low) & (source_5sels['sheared']['som_bin'][i] <= bin_high) & (source_5sels['sheared'][size_snr][i] <= np.median(source_5sels['sheared'][size_snr][i])) for i in range(5)]
        else:
                masks = [(source_5sels['sheared']['som_bin'][i] >= bin_low) & (source_5sels['sheared']['som_bin'][i] <= bin_high) & (source_5sels['sheared'][size_snr][i] > np.median(source_5sels['sheared'][size_snr][i])) for i in range(5)]
	
        source_bin = {}
        source_bin['ra'] = source['ra'][masks[0]]
        source_bin['dec'] = source['dec'][masks[0]]
        source_bin['e1'] = source['e1'][masks[0]]
        source_bin['e2'] = source['e2'][masks[0]]
        source_bin['psf_e1'] = source['psf_e1'][masks[0]]
        source_bin['psf_e2'] = source['psf_e2'][masks[0]]
        source_bin['snr'] = source['snr'][masks[0]]
        source_bin['size'] = source['size'][masks[0]]
        source_bin['bpz_mean'] = source['bpz_mean'][masks[0]]
        source_bin['bpz_zmc'] = source['bpz_zmc'][masks[0]]
        source_bin['Rgamma'] = source['Rgamma'][masks[0]]
	
        R11, _, _ = calibrator.calibrate('e_1', mask=masks)
        R22, _, _ = calibrator.calibrate('e_2', mask=masks)
        source_bin['Rmean'] = np.mean([R11, R22])
        source_bin['R11'] = calibrator.calibrate('e_1', return_full=True, mask=masks)[0]
        source_bin['R22'] = calibrator.calibrate('e_2', return_full=True, mask=masks)[0]
        print('Mean response redshift bin (%0.2f, %0.2f):'%(bin_low, bin_high), 
              source_bin['Rmean'], np.mean(source_bin['Rgamma']), np.mean(source_bin['R11']), np.mean(source_bin['R22']))
        return source_bin

    def run_responses_mean_notomo_size_snr(self, Rgamma, size_snr, cut, high_low, delta_gamma):
        """
        Computes responses when there is an extra selection on size or snr.
        For all the sources (no tomography).
        - Rgamma: Rgamma for the high or low selection in size or snr.
        - size_snr: string, either size or snr.
        - cut: median of size or snr for the whole source sample.
        - high_low: string, either high or low.
        - delta_gamma: value of the artificially applied shear to the images.      
        """

        e_ix = {}  # ellipticities for component i for a given selection s, divided by Delta gamma.
        # x: p, m
        components = ['1p', '1m', '2p', '2m']
        for comp in components:
            print(comp)
            e_ix_allbins = np.zeros(0)
            par_ix_allbins = np.zeros(0)  # size or snr in the selection ix, for instance, size_1p, size_2p, etc.
            for sbin in self.zbins['sbins']:
                # Appending responses source bin sbin.
                source_selection = pf.getdata(
                    self.paths['y1'] + 'metacal_sel_responses/metacal_sel_responses_sa%s_%s.fits' % (sbin[1], comp))
                print('dgamma is now dividing in the function that loads metacal and writes the file on it!')
                e_ix_allbins = np.append(e_ix_allbins, source_selection['Riisx'])
                par_ix_allbins = np.append(par_ix_allbins, source_selection['%s_ix' % size_snr])

            # For each component, mask the high or low part of e_ix
            if high_low == 'high':
                mask = np.where(par_ix_allbins > cut)
            if high_low == 'low':
                mask = np.where(par_ix_allbins <= cut)

            e_ix[comp] = np.mean(e_ix_allbins[mask])

        Rgamma_mean = np.mean(Rgamma)
        Rs_mean = self.compute_Rs(e_ix, delta_gamma)
        R_mean = self.save_responses_mean(Rgamma_mean, Rs_mean, 'notomo_%s_%s' % (size_snr, high_low))

        return R_mean

    def load_nzs(self, size_snr):
        """
	Loads redshift distributions for lenses, sources and source splits, low(l) and high(h).
	"""
        zl, nzl = np.loadtxt(self.paths['nz_lens'] + 'lens', unpack=True, usecols=(0, 1))
        zs, nzsl, nzsh = np.loadtxt(self.paths['nz_source_notomo_%s' % size_snr], unpack=True)
        nzsl = interpolate.interp1d(zs + self.source_nofz_pars['dzs', size_snr][0], nzsl, bounds_error=False,
                                    fill_value=0)(zs)
        nzsh = interpolate.interp1d(zs + self.source_nofz_pars['dzs', size_snr][1], nzsh, bounds_error=False,
                                    fill_value=0)(zs)

        return zl, nzl, zs, nzsl, nzsh

    def run(self, size_snr):
        lens_all, random_all, source_all, source_all_5sels, calibrator = self.load_data_or_sims()
        source = self.load_metacal_bin(source_all, source_all_5sels, calibrator, bin_low=self.zbins['s1'][0], bin_high=self.zbins['s4'][1])
        sourcel = self.load_metacal_bin_size_snr(source_all, source_all_5sels, calibrator, size_snr, 'low', bin_low=self.zbins['s1'][0], bin_high=self.zbins['s4'][1])
        Rl = sourcel['Rmean']

        sourceh = self.load_metacal_bin_size_snr(source_all, source_all_5sels, calibrator, size_snr, 'high', bin_low=self.zbins['s1'][0], bin_high=self.zbins['s4'][1])
        Rh = sourceh['Rmean']

        lens = lens_all[(lens_all['z'] > self.zbins['l1'][0]) & (lens_all['z'] < self.zbins['l1'][1])]
        random = random_all[(random_all['z'] > self.zbins['l1'][0]) & (random_all['z'] < self.zbins['l1'][1])]

        path_test = self.get_path_test(size_snr)
        make_directory(path_test)

        # Computing the measurements for the split halves, both around lenses and randoms
        theta, gtsl, gxsl, errsl, weightsl, npairsl = self.run_treecorr_jackknife(lens, sourcel, 'NG')
        self.save_runs(path_test, theta, gtsl, gxsl, errsl, weightsl, npairsl, False)
        gtlnum, gxlnum, wlnum = self.numerators_jackknife(gtsl, gxsl, weightsl)

        theta, gtsl_r, gxsl_r, errsl_r, weightsl_r, npairsl_r = self.run_treecorr_jackknife(random, sourcel, 'NG')
        self.save_runs(path_test, theta, gtsl_r, gxsl_r, errsl_r, weightsl_r, npairsl_r, True)
        gtlnum_r, gxlnum_r, wlnum_r = self.numerators_jackknife(gtsl_r, gxsl_r, weightsl_r)

        theta, gtsh, gxsh, errsh, weightsh, npairsh = self.run_treecorr_jackknife(lens, sourceh, 'NG')
        gthnum, gxhnum, whnum = self.numerators_jackknife(gtsh, gxsh, weightsh)

        theta, gtsh_r, gxsh_r, errsh_r, weightsh_r, npairsh_r = self.run_treecorr_jackknife(random, sourceh, 'NG')
        gthnum_r, gxhnum_r, whnum_r = self.numerators_jackknife(gtsh_r, gxsh_r, weightsh_r)

        # Combining measurements and responses to get gammat
        gtl_all = (gtlnum / wlnum) / Rl - (gtlnum_r / wlnum_r) / Rl
        gth_all = (gthnum / whnum) / Rh - (gthnum_r / whnum_r) / Rh

        bfl_all = self.compute_boost_factor(lens['jk'], random['jk'], wlnum, wlnum_r)
        gtl_all_boosted = bfl_all*(gtlnum / wlnum)/Rl -(gtlnum_r / wlnum_r)/Rl

        bfh_all = self.compute_boost_factor(lens['jk'], random['jk'], whnum, whnum_r)
        gth_all_boosted = bfh_all*(gthnum / whnum)/Rh -(gthnum_r / whnum_r)/Rh

        print('Measurements done. Saving n(z)s')

        np.savetxt(path_test + 'gtl_all', gtl_all)
        np.savetxt(path_test + 'gth_all', gth_all)
	
        self.process_run(gtl_all, theta, path_test, 'gtl')
        self.process_run(gth_all, theta, path_test, 'gth')
        self.process_run((gtlnum_r / wlnum_r) / Rl, theta, path_test, 'randomsl')
        self.process_run((gthnum_r / whnum_r) / Rh, theta, path_test, 'randomsh')
	
        self.process_run(bfl_all, theta, path_test, 'boost_factorl')
        self.process_run(bfh_all, theta, path_test, 'boost_factorh')
        self.process_run(gtl_all_boosted, theta, path_test, 'gtl_boosted')
        self.process_run(gth_all_boosted, theta, path_test, 'gth_boosted')
	
	
        #Computing the histograms and saving them
        if size_snr == 'size':
            hist_nofz_all_size = np.histogram(source['bpz_zmc'], normed=True, bins=200)
            hist_nofz_high_size = np.histogram(sourceh['bpz_zmc'], normed=True, bins=200)
            hist_nofz_low_size = np.histogram(sourcel['bpz_zmc'], normed=True, bins=200)
        else:
            hist_nofz_all_snr = np.histogram(source['bpz_zmc'], normed=True, bins=200)
            hist_nofz_high_snr = np.histogram(sourceh['bpz_zmc'], normed=True, bins=200)
            hist_nofz_low_snr = np.histogram(sourcel['bpz_zmc'], normed=True, bins=200)

        if 'combined_sample_fid' in self.config['lens_v']:
            nofz_lenses_witherr = np.random.normal(lens['z'],lens['zerr'])
        else:
            nofz_lenses_witherr = lens['zerr']

        hist_nofz_lenses_witherr = np.histogram(nofz_lenses_witherr, normed=True, bins=200)
        if size_snr == 'size':
            zstep_all_size=(hist_nofz_all_size[1][1]-hist_nofz_all_size[1][0])/2.
            zstep_high_size=(hist_nofz_high_size[1][1]-hist_nofz_high_size[1][0])/2.
            zstep_low_size=(hist_nofz_low_size[1][1]-hist_nofz_low_size[1][0])/2.
        else:
            zstep_all_snr=(hist_nofz_all_snr[1][1]-hist_nofz_all_snr[1][0])/2.
            zstep_high_snr=(hist_nofz_high_snr[1][1]-hist_nofz_high_snr[1][0])/2.
            zstep_low_snr=(hist_nofz_low_snr[1][1]-hist_nofz_low_snr[1][0])/2.

        zstep_lenses_witherr=(hist_nofz_lenses_witherr[1][1]-hist_nofz_lenses_witherr[1][0])/2.

        if size_snr == 'size':
            a=open(path_test + 'hist_n_of_z_all_size','w')
            for i in range(0,len(hist_nofz_all_size[1][0:-1])):
                a.write(str(hist_nofz_all_size[1][i]+zstep_all_size))
                a.write('\t')
                a.write(str(hist_nofz_all_size[0][i]))
                a.write('\n')
    
            a.close() 

            a=open(path_test + 'hist_n_of_z_high_size','w')
            for i in range(0,len(hist_nofz_high_size[1][0:-1])):
                a.write(str(hist_nofz_high_size[1][i]+zstep_high_size))
                a.write('\t')
                a.write(str(hist_nofz_high_size[0][i]))
                a.write('\n')

            a.close()

            a=open(path_test + 'hist_n_of_z_low_size','w')
            for i in range(0,len(hist_nofz_low_size[1][0:-1])):
                a.write(str(hist_nofz_low_size[1][i]+zstep_low_size))
                a.write('\t')
                a.write(str(hist_nofz_low_size[0][i]))
                a.write('\n')
    
            a.close()            

        else:
            a=open(path_test + 'hist_n_of_z_all_snr','w')
            for i in range(0,len(hist_nofz_all_snr[1][0:-1])):
                a.write(str(hist_nofz_all_snr[1][i]+zstep_all_snr))
                a.write('\t')
                a.write(str(hist_nofz_all_snr[0][i]))
                a.write('\n')
    
            a.close()

            a=open(path_test + 'hist_n_of_z_high_snr','w')
            for i in range(0,len(hist_nofz_high_snr[1][0:-1])):
                a.write(str(hist_nofz_high_snr[1][i]+zstep_high_snr))
                a.write('\t')
                a.write(str(hist_nofz_high_snr[0][i]))
                a.write('\n')

            a.close()

            a=open(path_test + 'hist_n_of_z_low_snr','w')
            for i in range(0,len(hist_nofz_low_snr[1][0:-1])):
                a.write(str(hist_nofz_low_snr[1][i]+zstep_low_snr))
                a.write('\t')
                a.write(str(hist_nofz_low_snr[0][i]))
                a.write('\n')

            a.close()

        a=open(path_test + 'hist_n_of_z_lenses_witherr','w')

        for i in range(0,len(hist_nofz_lenses_witherr[1][0:-1])):
            a.write(str(hist_nofz_lenses_witherr[1][i]+zstep_lenses_witherr))
            a.write('\t')
            a.write(str(hist_nofz_lenses_witherr[0][i]))
            a.write('\n')

        a.close()

        """
        np.savetxt(path_test + 'n_of_z_lenses', lens['z'])
            np.savetxt(path_test + 'n_of_z_lenses_err', lens['zerr'])
        if size_snr == 'size':
                np.savetxt(path_test + 'n_of_z_low_size', sourcel['bpz_zmc'])
                np.savetxt(path_test + 'n_of_z_high_size', sourceh['bpz_zmc'])
                np.savetxt(path_test + 'n_of_z_all_size', source['bpz_zmc'])
        else:
                np.savetxt(path_test + 'n_of_z_low_snr', sourcel['bpz_zmc'])
                np.savetxt(path_test + 'n_of_z_high_snr', sourceh['bpz_zmc'])
                np.savetxt(path_test + 'n_of_z_all_snr', source['bpz_zmc'])
        """
        print('All measurements done. Please, use the separate notebook to produce the plot. It still needs to be merged into the code.')
	
        """
        # Getting the data ratio using the simulations
        #sims, cov_sims = self.load_sims()
        
	import astropy.io.fits
	cov_theory=astropy.io.fits.open(self.paths['theory_size_all_covmat'])

	theory=[]
	for i in range(len(cov_theory[2].data)):
		theory.append(cov_theory[2].data[i][-2])

	theory=np.array(theory)

	ratio, err_ratio = self.ratio_from_sims(theta, gtl_all, gth_all, theory, cov_theory[1].data)

        # Load N(z)'s and corrects the mean using Cosmos calibration.
        #zl, nzl, zs, nzsl, nzsh = self.load_nzs(size_snr)
	zl, nzl = np.loadtxt(self.paths['hist_n_of_z_lenses_witherr_size'], unpack=True, usecols=(0, 1))
	zs, nzsl = np.loadtxt(self.paths['hist_n_of_z_low_size'], unpack=True, usecols=(0,1))
	zs, nzsh = np.loadtxt(self.paths['hist_n_of_z_high_size'], unpack=True, usecols=(0,1))

	nzsl = interpolate.interp1d(zs + self.source_nofz_pars['dzs', size_snr][0], nzsl, bounds_error=False,
                                    fill_value=0)(zs)
        nzsh = interpolate.interp1d(zs + self.source_nofz_pars['dzs', size_snr][1], nzsh, bounds_error=False,
                                    fill_value=0)(zs)

        # Computing inverse sigma_crit for the splits
        isch = functions.inv_sigma_crit_eff(zl, zs, nzl, nzsh)
        iscl = functions.inv_sigma_crit_eff(zl, zs, nzl, nzsl)

        # Computing the error on the theory prediction for the ratios, based on photo-z uncertainties
        # Shifting the N(z)'s for the splits up and down by the dzs_sigma parameter, and computing the theory ratio for each of these four cases
        shifts = [[-0.01, 0.], [0.01, 0.], [0., -0.01], [0., +0.01]]
        shifts = np.array(shifts) * self.source_nofz_pars['dzs_sigma']
        ratios = np.zeros(5)
        ratios[0] = isch / iscl
        sigmas = np.zeros(5)
        sigmas[0] = 0
        i = 1
        for s in shifts:
            nzsli = interpolate.interp1d(zs + s[0], nzsl, bounds_error=False, fill_value=0)
            nzsls = nzsli(zs)

            nzshi = interpolate.interp1d(zs + s[1], nzsh, bounds_error=False, fill_value=0)
            nzshs = nzshi(zs)

            ischs = functions.inv_sigma_crit_eff(zl, zs, nzl, nzshs)
            iscls = functions.inv_sigma_crit_eff(zl, zs, nzl, nzsls)
            ratios[i] = ischs / iscls
            sigmas[i] = abs(ratios[0] - ratios[i])
            i = i + 1

        # Taking the mean of the sigma up and down for each of the splits, and computing the total error by adding them in quadrature
        sigma_tot = np.sqrt(np.mean([sigmas[1:3]]) ** 2 + np.mean(sigmas[3:5]) ** 2)
        print ratio, err_ratio
        print '%0.2f +- %0.2f' % (ratios[0], sigma_tot)

        # Saving the data ratio
        result_data = [ratio, err_ratio]
        result_theory = [ratios[0], sigma_tot]
        self.save_size_snr(path_test, size_snr, result_data, result_theory)

        """
    
    def plot(self):
        """
        plt.rc('text', usetex=self.plotting['latex'])
        plt.rc('font', family='serif')

        cmap = self.plotting['cmap']
        cmap_step = 0.25
        lss = ['-', ':']
        labels_c = [self.plotting['catname']]
        tests = ['snr', 'size']
        labels_t = [r'S/N', r'Size']

        # Text
        titles_l = r'$(0.15 < z_l < 0.30)$'
        titles_s = r'$0.20 < z_s < 1.30$'

        fig, ax = plt.subplots(2, 2, figsize=(8, 6), sharey='row', sharex=False, gridspec_kw={'height_ratios': [1, 1]})
        fig.subplots_adjust(hspace=0.32, wspace=0.0)

        for t, size_snr in enumerate(tests):
            color = plt.get_cmap(cmap)(0)
            # Load N(z)'s and corrects the mean using Cosmos calibration.
            zl, nzl, zs, nzsl, nzsh = self.load_nzs(size_snr)
            ax[0][t].fill_between(zl, 0, nzl / np.max(nzl) * 2.04696244339, color='gray', alpha=0.2)
            ax[0][t].set_xlim(0, 1.8)
            ax[0][t].set_xlabel('Redshift')
            ax[0][t].plot(zs + self.source_nofz_pars['dzs', size_snr][0], nzsl, ls='-', lw=1.5, c=color,
                          label='Low %s' % (labels_c[0]))
            ax[0][t].plot(zs + self.source_nofz_pars['dzs', size_snr][1], nzsh, ls='-.', lw=1.5, c=color,
                          label='High %s' % (labels_c[0]))
            ax[0][t].set_title('%s split' % labels_t[t])

        ax[0][0].plot([], [], color='gray', alpha=0.2, linewidth=10, label='redMaGiC ' + titles_l)
        ax[0][0].set_ylabel('$n(z)$', fontsize=14)
        c1 = plt.get_cmap(cmap)(cmap_step * 1.5)
        c2 = plt.get_cmap(cmap)(cmap_step * 3)
        c3 = 'k'

        handles, labels = ax[0][0].get_legend_handles_labels()

        # Sort the legend by the labels names
        # ax[0][1].legend(handles, labels, frameon=False,loc='best', prop={'size':9})
        # ax[0][0].set_xlim(0,1.79)
        ax[0][0].set_xlim(0, 1.79)
        ax[0][0].set_ylim(0, 3.)
        ax[0][1].legend(frameon=False, loc='best', prop={'size': 9})

        def compute_sigmas(q1, err1, q2, err2):
            # q1, err1: quantity and corresponding error
            print 'q1,err1,q2,err2:', q1, err1, q2, err2
            return np.abs(q1 - q2) / np.sqrt(err1 ** 2 + err2 ** 2)

        for t, size_snr in enumerate(tests):
            path_test = self.get_path_test(size_snr)
            ratio, err_ratio, invsc, err_invsc = self.load_size_snr(path_test, size_snr)
            sigma = compute_sigmas(ratio, err_ratio, invsc, err_invsc)

            color, x = c1, 1
            ax[1][t].errorbar(x, ratio, err_ratio, fmt='o', color=color, mec=color, label='%s' % (labels_c[0]))
            x_lin = [x - 0.2, x, x + 0.2]
            ax[1][t].fill_between(x_lin, np.array([invsc - err_invsc for i in x_lin]),
                                  np.array([invsc + err_invsc for i in x_lin]), alpha=0.4, edgecolor=color,
                                  facecolor=color,
                                  label='$\Sigma_{\mathrm{crit,eff}}^{-1,\,\mathrm{high}}/ \Sigma_{\mathrm{crit,eff}}^{-1,\,\mathrm{low}}$ %s' % (
                                      labels_c[0]))

            ax[1][t].set_xlim(0.5, 1.5)
            ax[1][t].set_ylim(0.4, 1.6)
            ax[1][t].set_xticklabels([])
            ax[1][t].set_xticks([])

        ax[1][0].legend(frameon=False, prop={'size': 10.5}, loc='best', numpoints=1)
        ax[1][0].set_ylabel(r'$\gamma^{\mathrm{high}}_t /\ \gamma^{\mathrm{low}}_t $', fontsize=14)

        self.save_plot('snr_size')
        """

class init_ggl_class(GGL):

    def __init__(self, basic, config, paths, zbins, downsample=1):
        GGL.__init__(self, basic, config, paths)
        self.zbins = zbins

        if self.basic['mode'] == 'data':
            lens_all, random_all, source_all, source_all_5sels, calibrator = self.load_data_or_sims()

        for sbin in self.zbins['sbins']:
            print('Running measurement for source %s.' % sbin)

        if self.basic['mode'] == 'data':
            self.source_all = self.load_metacal_bin(source_all, source_all_5sels, calibrator,
                                                    zlim_low=self.zbins['sys'][0], zlim_high=self.zbins['sys'][1])
            self.R = self.source_all['Rmean']
            print(np.min(self.source_all['bpz_mean'])), (np.max(self.source_all['bpz_mean']))

        if downsample > 1:
            ind_ds = np.random.randint(0, len(lens_all['z']), int(len(lens_all['z']) / downsample))
            lens_all = lens_all[ind_ds]

        self.lens = lens_all

        if downsample > 1:
            ind_ds = np.random.randint(0, len(random_all['z']), int(len(random_all['z']) / downsample))
            random_all = random_all[ind_ds]

        self.random = random_all

        keys_new_dict = ['ra', 'dec', 'e1', 'e2', 'bpz_zmc']
        self.source = {}

        for keys in keys_new_dict:
            if downsample > 1:
                ind_ds = np.random.randint(0, len(self.source_all[keys]), int(len(self.source_all[keys]) / downsample))
                self.source[keys] = self.source_all[keys][ind_ds]
            else:
                self.source[keys] = self.source_all[keys]
        del source_all, random_all, lens_all


class TestSysMaps():
    """
    SubClass to test if the tangential shear signal has no dependence on observational conditions such as seeing, airmass etc, in each band griz.
    Using the first lens bin and all sources, as in Y1.
    The variable map can be: 'airmass', 'fwhm', 'maglimit', 'skybrite'. We iterate over them.
    The variable band can be: 'g', 'r', 'i', 'z'. We iterate over them. In Y1 we only used r band in the end, because it was used by im3shape.
    """

    def __init__(self, ggl_class, config, paths, zbins, plotting, source_nofz_pars, sysmaps):
        self.zbins = zbins
        self.plotting = plotting
        self.source_nofz_pars = source_nofz_pars
        self.sysmaps = sysmaps
        self.ggl = ggl_class
        self.paths = paths

    def get_path_test(self, map, band):
        return os.path.join(self.paths['runs_config'], 'systematics_maps', map, band) + '/'

    def save_systematics_maps_ratios(self, path_test, result_data, result_theory):
        save_suffix = ''
        if self.sysmaps['separate_jk']:
            njk = self.sysmaps['separate_jk_njk']
            save_suffix = save_suffix + '_separate_jk_njk_' + str(njk)

        np.savetxt(path_test + 'data' + save_suffix, result_data, header='ratio err_ratio')
        np.savetxt(path_test + 'theory' + save_suffix, result_theory, header='ratio err_ratio')

    def load_systematics_maps_ratios(self, path_test):
        save_suffix = ''
        if self.sysmaps['separate_jk']:
            njk = self.sysmaps['separate_jk_njk']
            save_suffix = save_suffix + '_separate_jk_njk_' + str(njk)

        data, data_err = np.loadtxt(path_test + 'data' + save_suffix, unpack=True)
        theory, theory_err = np.loadtxt(path_test + 'theory' + save_suffix, unpack=True)
        return data, data_err, theory, theory_err


    def load_theory(self):
        dataf = pf.open(
            '/global/project/projectdirs/des/shivamp/ggl_results/data_cov_sys/gammat_theory_FLASKcosmology_covmatG_upsample10_halfarea.fits')
        data = dataf[2].data['VALUE']
        cov = dataf[1].data
        return data, cov

    def ratio_from_theory(self, theta, gtl_all, gth_all, sims, cov_sims):

        ratioA_all = np.zeros(len(gth_all))
        chi2fit_h = np.zeros(len(gth_all))
        chi2fit_l = np.zeros(len(gth_all))

        mask = (theta >= source_nofz_pars['thetamin'])
        for i in range(len(gth_all)):
            Ah, chi2fit_h[i], _ = functions.minimize_chi2_fit_amplitude((cov_sims[mask].T)[mask].T, gth_all[i][mask],
                                                                        sims[mask])
            Al, chi2fit_l[i], _ = functions.minimize_chi2_fit_amplitude((cov_sims[mask].T)[mask].T, gtl_all[i][mask],
                                                                        sims[mask])
            ratioA_all[i] = Ah / Al

        ratioA_mean = np.mean(ratioA_all, axis=0)
        COV_A = (len(ratioA_all) - 1) * np.cov(ratioA_all, bias=True)
        err_A = np.sqrt(COV_A)
        return ratioA_mean, err_A


    def save_gammat(self, all, theta, path_test, end):
        """
        From the jackknife measurements in all jackknife regions but all, constructs covariance, mean and stats.
        Saves them into file.
        all: gt_all or gx_all.
        theta: in arcmin.
        path_test: where to save the files.
        end: string to save the files: gt, gx, randoms, etc.
        """
        mean = np.mean(all, axis=0)
        cov = functions.covariance(all, mean)
        err = np.sqrt(np.diag(cov))
        chi2 = np.dot(mean.T, np.dot(np.linalg.inv(cov), mean))
        ndf = len(mean)

        std = np.sqrt((len(all) - 1.)) * np.std(all, axis=0)

        # Hartlap factor
        N = len(all)  # number of jackknife regions
        p = len(mean)  # number of angular bins
        factor = (N - p - 2) / float(N - 1)
        chi2_hartlap = chi2 * factor

        stats = np.array([chi2_hartlap, ndf])

        np.savetxt(path_test + 'mean_%s' % end, zip(theta, mean, err), header='th, %s, err_%s' % (end, end))
        np.savetxt(path_test + 'cov_%s' % end, cov)
        np.savetxt(path_test + 'all_%s' % end, all,
                   header='%s (sum of %s from all patches except for one, different each time)' % (end, end))
        np.savetxt(path_test + 'null_chi2_%s' % end, stats.reshape(1, stats.shape[0]),
                   fmt='%0.1f  %d', header='chi2_hartlap  ndf')

        return theta, mean, err

    def visualize_map(self, pix, signal, map, band, nside, nested_bool, name):
        """
        Plots the map and saves it.
        """
        ma = np.array([-1.6375 * 10 ** 30] * hp.nside2npix(nside))
        naturals = np.arange(0, len(ma))
        pix = np.in1d(naturals, pix)
        ma[pix] = signal
        # hp.gnomview(ma, min= mini, max=maxi, rot =  (73, -52, 0), xsize = 750 , title = '%s i'%map ,  notext = True, cmap = self.plotting['cmap'], nest = nested_bool)
        hp.mollview(ma, title='%s %s' % (map, band), notext=True, cmap=self.plotting['cmap'], nest=nested_bool)
        path = self.paths['plots'] + 'systematics_maps/'
        make_directory(path)
        plt.savefig(path + '%s_%s%s.pdf' % (map, band, name))

    def load_systematics_map_y3(self, map, band, nside, nested_bool):
        '''
        Loads the systematics map, splits into high and low parts, and plots all maps.
        nested_bool: True if nest, False if ring.
        Returns: pixels corresponding the low half and the high half, for each map.
        '''
        path = (self.paths['y3_sysmap'])
        sys_map = pf.open(path + 'y3a2_%s' % band + '_o.%s' % nside + '_t.32768_' + self.sysmaps[map] + '.fits')[1].data
        pix = sys_map['PIXEL']
        sig = sys_map['SIGNAL']
        self.visualize_map(pix, sig, map, band, nside, nested_bool, '')

        mask_low = sig <= np.median(sig)
        mask_hi = sig > np.median(sig)

        pix_low = pix[mask_low]
        sig_low = sig[mask_low]

        pix_hi = pix[mask_hi]
        sig_hi = sig[mask_hi]

        self.visualize_map(pix_low, sig_low, map, band, nside, nested_bool, '_low')
        self.visualize_map(pix_hi, sig_hi, map, band, nside, nested_bool, '_high')

        return pix_low, pix_hi

    def load_normedsystematics_map_y3(self, map, band, nside, nested_bool):
        '''
        Loads the systematics map, splits into high and low parts, and plots all maps.
        nested_bool: True if nest, False if ring.
        Returns: pixels corresponding the low half and the high half, for each map.
        '''
        path = (self.paths['y3_sysmap'])
        sys_map = pf.open(path + 'y3a2_%s' % band + '_o.%s' % nside + '_t.32768_' + self.sysmaps[map] + '.fits')[1].data
        pix = sys_map['PIXEL']
        sig = sys_map['SIGNAL']
        self.visualize_map(pix, sig, map, band, nside, nested_bool, '')

        sig_normed = sig / np.mean(sig)

        return pix, sig_normed

    def load_nzs_y3(self, bpz_zmc, zmin, zmax, mask_low, mask_high):

        nzlf = np.loadtxt('/global/u1/s/spandey/xcorr/notebooks/nz_lens_y3_bin1.txt')
        zl, nzl = nzlf[:, 0], nzlf[:, 1]

        nzbins_total = 150

        delta_z = (zmax - zmin) / nzbins_total
        zarray_all = np.linspace(zmin, zmax, nzbins_total)
        zarray_edges = (zarray_all[1:] + zarray_all[:-1]) / 2.
        zarray = zarray_all[1:-1]

        hist_zmc_low, bin_edges = np.histogram(bpz_zmc[mask_low], bins=zarray_edges)

        hist_zmc_high, bin_edges = np.histogram(bpz_zmc[mask_high], bins=zarray_edges)

        nzsl = hist_zmc_low / (np.sum(hist_zmc_low) * delta_z)

        nzsh = hist_zmc_high / (np.sum(hist_zmc_high) * delta_z)

        return zl, nzl, zarray, nzsl, nzsh

    def radec_to_thetaphi(self, ra, dec):
        """
        Converts ra and dec in degrees to theta and phi.
        Returns theta and phi in radians.
        """
        theta = (90. - dec) * np.pi / 180.
        phi = ra * np.pi / 180.
        return theta, phi

    def plot_uw_ww(self, uw_array, w_array, save_dir, save_suffix):
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        theta_uw, mean_uw, err_uw = uw_array[0], uw_array[1], uw_array[2]
        theta_w, mean_w, err_w = w_array[0], w_array[1], w_array[2]

        ax.errorbar(theta_uw, theta_uw * mean_uw, theta_uw * err_uw, label='Unweighted')
        ax.errorbar(theta_w * 1.02, theta_w * mean_w, theta_w * err_w, label='Weighted')

        ax.set_xscale('log')
        ax.set_yscale('log')

        ax.set_xlabel(r'$\theta$', size=20)
        ax.set_ylabel(r'$\theta \times \gamma_t$', size=20)

        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.tick_params(axis='both', which='minor', labelsize=15)

        ax.legend(fontsize=20, frameon=False)

        plt.savefig(save_dir + 'uw_ww_gammat_' + save_suffix + 'png')

    def run_wwdd(self, maps, bands):
        """
            Runs weighted gglensing measurment for all maps and bands.
        """

        for map in maps:
            print('Running map %s...' % map)
            for band in bands:
                print('Band %s' % band)

                path_test = self.get_path_test(map, band)
                make_directory(path_test)

                # Load and split the systematics map
                pix_all_sys, sys_map_normed = self.load_normedsystematics_map_y3(map, band, self.sysmaps['nside'],
                                                                                 self.sysmaps['nested_bool'])

                # Building lenses masks, low and high
                theta_l, phi_l = self.radec_to_thetaphi(self.ggl.lens['ra'], self.ggl.lens['dec'])
                pix_all_l = hp.pixelfunc.ang2pix(self.sysmaps['nside'], theta_l, phi_l, self.sysmaps['nested_bool'])
                mask_sys = np.in1d(pix_all_l, pix_all_sys)
                self.lens_masked = self.lens[mask_sys]
                theta_l, phi_l = self.radec_to_thetaphi(self.lens_masked['ra'], self.lens_masked['dec'])
                pix_all_l = hp.pixelfunc.ang2pix(self.sysmaps['nside'], theta_l, phi_l, self.sysmaps['nested_bool'])
                ind_weight = np.searchsorted(pix_all_sys, pix_all_l)
                weight_l = sys_map_normed[ind_weight]
                self.lens_weighted = copy.deepcopy(self.lens_masked)
                self.lens_weighted['w'] = weight_l

                # Building randoms masks, low and high
                theta_r, phi_r = self.radec_to_thetaphi(self.ggl.random['ra'], self.ggl.random['dec'])
                pix_all_r = hp.pixelfunc.ang2pix(self.sysmaps['nside'], theta_r, phi_r, self.sysmaps['nested_bool'])
                mask_sys = np.in1d(pix_all_r, pix_all_sys)
                self.random_masked = self.random[mask_sys]

                # Building sources masks, low and high
                theta_s, phi_s = self.radec_to_thetaphi(self.ggl.source['ra'], self.ggl.source['dec'])
                pix_all_s = hp.pixelfunc.ang2pix(self.sysmaps['nside'], theta_s, phi_s, self.sysmaps['nested_bool'])
                mask_sys = np.in1d(pix_all_s, pix_all_sys)

                keys_new_dict = ['ra', 'dec', 'e1', 'e2']
                self.source_masked = {}

                for keys in keys_new_dict:
                    self.source_masked[keys] = self.source[keys][mask_sys]

                ipdb.set_trace()
                # Computing the measurements for the split halves, both around lenses and randoms
                print('Lenses, unweighted')
                self.theta, self.gtsl, self.gxsl, self.errsl, self.weightsl, self.npairsl = self.ggl.run_treecorr_jackknife(
                    self.lens_masked, self.source_masked, 'NG')
                self.gtlnum, self.gxlnum, self.wlnum = self.ggl.numerators_jackknife(self.gtsl, self.gxsl,
                                                                                     self.weightsl)

                print('Randoms, low.')
                self.theta, self.gtsl_r, self.gxsl_r, self.errsl_r, self.weightsl_r, self.npairsl_r = self.ggl.run_treecorr_jackknife(
                    self.random_masked,
                    self.source_masked, 'NG')
                self.gtlnum_r, self.gxlnum_r, self.wlnum_r = self.ggl.numerators_jackknife(self.gtsl_r, self.gxsl_r,
                                                                                           self.weightsl_r)

                print('Lenses, weighted')
                self.theta, self.gtsh, self.gxsh, self.errsh, self.weightsh, self.npairsh = self.ggl.run_treecorr_jackknife(
                    self.lens_weighted, self.source_masked, 'NG')
                self.gthnum, self.gxhnum, self.whnum = self.ggl.numerators_jackknife(self.gtsh, self.gxsh,
                                                                                     self.weightsh)

                # Combining measurements and responses to get gammat
                gtl_all = (self.gtlnum / self.wlnum) / self.R - (self.gtlnum_r / self.wlnum_r) / self.R
                gth_all = (self.gthnum / self.whnum) / self.R - (self.gtlnum_r / self.wlnum_r) / self.R

                save_suffix = ''
                if self.sysmaps['separate_jk']:
                    njk = self.sysmaps['separate_jk_njk']
                    save_suffix = save_suffix + '_separate_jk_njk_' + str(njk)

                np.savez(
                    '/global/u1/s/spandey/ggl_results/gt_weighted_unweighted_npairs_zl_0.15_0.35_zs_0.2_1.2_band_' + band + '_' + str(
                        map) + save_suffix + '.npz',
                    gt_weighted=gth_all, gt_unweighted=gtl_all, theta=theta,
                    npairsl_unweighted=npairsl, npairs_rand=npairsl_r, npairs_weighted=npairsh)

                theta_uw, mean_uw, err_uw = self.save_gammat(gtl_all, theta, '/global/u1/s/spandey/ggl_results/',
                                                             'unweighted_band_' + band + '_' + str(map))
                theta_w, mean_w, err_w = self.save_gammat(gth_all, theta, '/global/u1/s/spandey/ggl_results/',
                                                          'weighted_band_' + band + '_' + str(map))

                self.plot_uw_ww([theta_uw, mean_uw, err_uw], [theta_w, mean_w, err_w],
                                '/global/u1/s/spandey/xcorr/plots/weighted_unweighted/',
                                'band_' + band + '_' + str(map))

    def run_y3(self, maps, bands):

        R = self.ggl.R
        for map in maps:
            print('Running map %s...' % map)
            for band in bands:
                print('Band %s' % band)

                path_test = self.get_path_test(map, band)
                make_directory(path_test)

                # Load and split the systematics map
                pix_low, pix_hi = self.load_systematics_map_y3(map, band, self.sysmaps['nside'],
                                                               self.sysmaps['nested_bool'])

                # Building lenses masks, low and high
                theta_l, phi_l = self.radec_to_thetaphi(self.ggl.lens['ra'], self.ggl.lens['dec'])
                pix_all_l = hp.pixelfunc.ang2pix(self.sysmaps['nside'], theta_l, phi_l, self.sysmaps['nested_bool'])
                maskl_low = np.in1d(pix_all_l, pix_low)
                maskl_hi = np.in1d(pix_all_l, pix_hi)

                # Building randoms masks, low and high
                theta_r, phi_r = self.radec_to_thetaphi(self.ggl.random['ra'], self.ggl.random['dec'])
                pix_all_r = hp.pixelfunc.ang2pix(self.sysmaps['nside'], theta_r, phi_r, self.sysmaps['nested_bool'])
                maskr_low = np.in1d(pix_all_r, pix_low)
                maskr_hi = np.in1d(pix_all_r, pix_hi)

                # Building sources masks, low and high
                theta_s, phi_s = self.radec_to_thetaphi(self.ggl.source['ra'], self.ggl.source['dec'])
                pix_all_s = hp.pixelfunc.ang2pix(self.sysmaps['nside'], theta_s, phi_s, self.sysmaps['nested_bool'])
                masks_low = np.in1d(pix_all_s, pix_low)
                masks_hi = np.in1d(pix_all_s, pix_hi)

                # ipdb.set_trace()

                if self.sysmaps['separate_jk']:
                    njk = self.sysmaps['separate_jk_njk']
                    lens_ra_l, lens_dec_l = self.ggl.lens['ra'][maskl_low], self.ggl.lens['dec'][maskl_low]
                    jkobj_map_l = kmeans_radec.kmeans_sample(np.transpose([lens_ra_l, lens_dec_l]), njk, maxiter=200)
                    lens_low_jk = jkobj_map_l.find_nearest(np.transpose([lens_ra_l, lens_dec_l]))
                    lens_low = self.ggl.lens[maskl_low]
                    lens_low['jk'] = lens_low_jk

                    lens_ra_h, lens_dec_h = self.ggl.lens['ra'][maskl_hi], self.ggl.lens['dec'][maskl_hi]
                    jkobj_map_h = kmeans_radec.kmeans_sample(np.transpose([lens_ra_h, lens_dec_h]), njk, maxiter=200)
                    lens_hi_jk = jkobj_map_h.find_nearest(np.transpose([lens_ra_h, lens_dec_h]))
                    lens_hi = self.ggl.lens[maskl_hi]
                    lens_hi['jk'] = lens_hi_jk

                    random_low = self.ggl.random[maskr_low]
                    random_low_jk = jkobj_map_l.find_nearest(np.transpose([random_low['ra'], random_low['dec']]))
                    random_low['jk'] = random_low_jk

                    random_hi = self.ggl.random[maskr_hi]
                    random_hi_jk = jkobj_map_h.find_nearest(np.transpose([random_hi['ra'], random_hi['dec']]))
                    random_hi['jk'] = random_hi_jk
                else:
                    lens_low = self.ggl.lens[maskl_low]
                    lens_hi = self.ggl.lens[maskl_hi]
                    random_low = self.ggl.random[maskr_low]
                    random_hi = self.ggl.random[maskr_hi]

                # Computing the measurements for the split halves, both around lenses and randoms
                print('Lenses, low.')
                theta, gtsl, gxsl, errsl, weightsl, npairsl = self.ggl.run_treecorr_jackknife(lens_low,
                                                                                              self.ggl.source, 'NG')
                gtlnum, gxlnum, wlnum = self.ggl.numerators_jackknife(gtsl, gxsl, weightsl)

                print('Randoms, low.')
                theta, gtsl_r, gxsl_r, errsl_r, weightsl_r, npairsl_r = self.ggl.run_treecorr_jackknife(random_low,
                                                                                                        self.ggl.source,
                                                                                                        'NG')
                gtlnum_r, gxlnum_r, wlnum_r = self.ggl.numerators_jackknife(gtsl_r, gxsl_r, weightsl_r)

                print('Lenses, high.')
                theta, gtsh, gxsh, errsh, weightsh, npairsh = self.ggl.run_treecorr_jackknife(lens_hi,
                                                                                              self.ggl.source, 'NG')
                gthnum, gxhnum, whnum = self.ggl.numerators_jackknife(gtsh, gxsh, weightsh)

                print('Randoms, high.')
                theta, gtsh_r, gxsh_r, errsh_r, weightsh_r, npairsh_r = self.ggl.run_treecorr_jackknife(random_hi,
                                                                                                        self.ggl.source,
                                                                                                        'NG')
                gthnum_r, gxhnum_r, whnum_r = self.ggl.numerators_jackknife(gtsh_r, gxsh_r, weightsh_r)

                # Combining measurements and responses to get gammat
                gtl_all = (gtlnum / wlnum) / R - (gtlnum_r / wlnum_r) / R
                gth_all = (gthnum / whnum) / R - (gthnum_r / whnum_r) / R

                save_suffix = ''
                if self.sysmaps['separate_jk']:
                    njk = self.sysmaps['separate_jk_njk']
                    save_suffix = save_suffix + '_separate_jk_njk_' + str(njk)

                np.savez('/global/u1/s/spandey/ggl_results/gtl_gth_npairs_zl_0.15_0.35_zs_0.2_1.2_band_' + band + '_' + \
                         str(map) + save_suffix + '.npz', gth=gth_all, gtl=gtl_all, theta=theta,
                         npairsl=npairsl, npairsl_r=npairsl_r, npairsh=npairsh, npairsh_r=npairsh_r)

                sims, cov_sims = self.load_theory()
                ratio, err_ratio = self.ratio_from_theory(theta, gtl_all, gth_all, sims, cov_sims)

                # Load N(z)'s
                zl, nzl, zs, nzsl, nzsh = self.load_nzs_y3(self.ggl.source['bpz_zmc'], zbins['sys'][0], zbins['sys'][1],
                                                           masks_low, masks_hi)

                # Computing inverse sigma_crit for the splits
                isch = functions.inv_sigma_crit_eff(zl, zs, nzl, nzsh)
                iscl = functions.inv_sigma_crit_eff(zl, zs, nzl, nzsl)

                # Computing the error on the theory prediction for the ratios, based on photo-z uncertainties
                # Shifting the N(z)'s for the splits up and down by the dzs_sigma parameter, and computing the theory ratio for each of these four cases
                shifts = [[-0.01, 0.], [0.01, 0.], [0., -0.01], [0., +0.01]]
                shifts = np.array(shifts) * self.source_nofz_pars['dzs_sigma']
                ratios = np.zeros(5)
                ratios[0] = isch / iscl
                sigmas = np.zeros(5)
                sigmas[0] = 0
                i = 1
                for s in shifts:
                    nzsli = interpolate.interp1d(zs + s[0], nzsl, bounds_error=False, fill_value=0)
                    nzsls = nzsli(zs)

                    nzshi = interpolate.interp1d(zs + s[1], nzsh, bounds_error=False, fill_value=0)
                    nzshs = nzshi(zs)

                    ischs = functions.inv_sigma_crit_eff(zl, zs, nzl, nzshs)
                    iscls = functions.inv_sigma_crit_eff(zl, zs, nzl, nzsls)
                    ratios[i] = ischs / iscls
                    sigmas[i] = abs(ratios[0] - ratios[i])
                    i = i + 1

                # Taking the mean of the sigma up and down for each of the splits, and computing the total error by adding them in quadrature
                sigma_tot = np.sqrt(np.mean([sigmas[1:3]]) ** 2 + np.mean(sigmas[3:5]) ** 2)
                print(ratio, err_ratio)
                print('%0.2f +- %0.2f' % (ratios[0], sigma_tot))

                # Saving the data ratio
                result_data = [ratio, err_ratio]
                result_theory = [ratios[0], sigma_tot]
                self.save_systematics_maps_ratios(path_test, result_data, result_theory)

    def plot_y3(self, bands):

        labels_cshort = 'Metacal'
        fontsize = 16
        maps = ['maglimit', 'skybrite', 'airmass', 'fwhm']
        #         maps = ['maglimit',  'skybrite']
        # maps = ['airmass', 'fwhm']
        # bands = 'r'
        cmap = self.plotting['cmap']
        cmap_step = 0.25
        c1 = plt.get_cmap(cmap)(cmap_step * 1.5)
        c2 = plt.get_cmap(cmap)(cmap_step * 3)

        fig, ax = plt.subplots(2, 2, figsize=(6, 6), sharey=True, sharex=True)
        fig.subplots_adjust(hspace=0.0, wspace=0.00)

        for m, map in enumerate(maps):
            if m % 2 == 0: c = 0
            if m % 2 == 1: c = 1
            if m == 0 or m == 1: f = 0
            if m == 2 or m == 3: f = 1

            x = 1
            j = 0
            for band in bands:
                color = plt.get_cmap(cmap)(cmap_step * 1 * (j + 1))
                path_test = self.get_path_test(map, band)
                data, data_err, theory, theory_err = self.load_systematics_maps_ratios(path_test)

                ax[f][c].errorbar(x - 0.1 + 0.1 * j + 1 * j, data, data_err, fmt='o', color=color, mec=color,
                                  label='band ' + band)

                x_lin = [x - 0.2 + 1 * j, x + 1 * j, x + 0.2 + 1 * j]

                if band == 'r':
                    ax[f][c].fill_between(x_lin, np.array([theory - theory_err for i in x_lin]),
                                          np.array([theory + theory_err for i in x_lin]), alpha=0.4, edgecolor=color,
                                          facecolor=color,
                                          label=r'$\Sigma_{\mathrm{crit,eff}}^{-1,\,\mathrm{high}}/ \Sigma_{\mathrm{crit,eff}}^{-1,\,\mathrm{low}}$')
                else:
                    ax[f][c].fill_between(x_lin, np.array([theory - theory_err for i in x_lin]),
                                          np.array([theory + theory_err for i in x_lin]), alpha=0.4, edgecolor=color,
                                          facecolor=color)
                j += 1

            ax[f][c].set_xlim(0.0, x + 1 * j + 1)
            ax[f][c].set_xticklabels([])
            ax[f][c].set_xticks([])
            ax[f][c].set_ylim(0.5, 1.5)
            ax[f][0].set_ylabel(r'$\gamma^{\mathrm{high}}_t /\ \gamma^{\mathrm{low}}_t $', fontsize=fontsize)
            ax[f][c].text(0.5, 0.9, map,
                          horizontalalignment='center', verticalalignment='center', transform=ax[f][c].transAxes,
                          fontsize=12, color='k')

        ax[1][0].legend(frameon=False, loc='best', numpoints=1, fontsize=9.8)

        save_suffix = ''
        if self.sysmaps['separate_jk']:
            njk = self.sysmaps['separate_jk_njk']
            save_suffix = save_suffix + '_separate_jk_njk_' + str(njk)

        self.ggl.save_plot('systematics_maps' + ''.join(bands) + save_suffix)
