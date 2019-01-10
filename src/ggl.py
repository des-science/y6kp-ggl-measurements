import psutil
import os
import ipdb
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
import itertools
import signal
import twopoint
from scipy import interpolate
import functions
from info import blind, paths, config, zbins, plotting, source_nofz_pars, sysmaps, mode, filename_mastercat, plot_blinded
import sys
import yaml
import subprocess
sys.path.append('../../destest/')
import destest

def make_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


plt.rc('text', usetex=False)
plt.rc('font', family='serif')


class GGL(object):
    """
    Basic class that has all the functions shared for several tests.
    """

    def __init__(self, config, paths):
        self.config = config
        self.paths = paths

    def load_metacal(self):
        """
        Loads metacal data for Y3 catalog using h5 interface.
        Combines with Gold and BPZ.
        Returns: dictionary for the sources, with all the relevant columns.
        """

        mcal_file = self.paths['yaml'] + 'destest_mcal.yaml'
        params_mcal = yaml.load(open(mcal_file))
        params_mcal['param_file'] = mcal_file
        params_mcal['filename'] = filename_mastercat
        source_mcal = destest.H5Source(params_mcal)
        source_selector = destest.Selector(params_mcal,source_mcal)
        source_calibrator = destest.MetaCalib(params_mcal,source_selector)

        gold_file = self.paths['yaml'] + 'destest_gold.yaml'
        params_gold = yaml.load(open(gold_file))
        params_gold['param_file'] = gold_file
        params_gold['filename'] = filename_mastercat
        source_gold = destest.H5Source(params_gold)
        gold_selector = destest.Selector(params_gold,source_gold,inherit=source_selector)

        param_file = self.paths['yaml'] + './destest_pz.yaml'
        params_pz = yaml.load(open(param_file))
        params_pz['filename'] = filename_mastercat
        source_pz = destest.H5Source(params_pz)
        pz_selector = destest.Selector(params_pz, source_pz, inherit=source_selector)

        source = {}
        source['ra'] = gold_selector.get_col('ra')[0]
        source['dec'] = gold_selector.get_col('dec')[0]
        if 'v1' in self.config['mastercat_v']:
            source['e1'] = source_selector.get_col('e1')[0]
            source['e2'] = source_selector.get_col('e2')[0]
        if 'v2' in self.config['mastercat_v']:
            source['e1'] = source_selector.get_col('e_1')[0]
            source['e2'] = source_selector.get_col('e_2')[0]
        source['psf_e1'] = source_selector.get_col('psf_e1')[0]
        source['psf_e2'] = source_selector.get_col('psf_e2')[0]
        source['snr'] = source_selector.get_col('snr')[0]
        source['size'] = source_selector.get_col('T')[0]
        source['bpz_mean'] = pz_selector.get_col('zmean_sof')
        source['bpz_zmc'] = pz_selector.get_col('zmc_sof')

        calibrator = destest.MetaCalib(params_mcal, source_selector)
        if 'v1' in self.config['mastercat_v']:
            R11, _, _ = calibrator.calibrate('e1')
            R22, _, _ = calibrator.calibrate('e2')
        if 'v2' in self.config['mastercat_v']:
            R11, _, _ = calibrator.calibrate('e_1')
            R22, _, _ = calibrator.calibrate('e_2')

        source['Rmean'] = np.mean([R11, R22])
        print 'Response full sample', source['Rmean']
        return source, calibrator

    def load_metacal_bin(self, source, calibrator, zlim_low, zlim_high):
        """
        source: dictionary containing relevant columns for the sources, with the baseline selection applied already.
        calibrator: class to compute the response. Taken from baseline selection.
        zlim_low, zlim_high: limits to select the tomographic bin.
        Obtains 5 masks (unsheared, sheared 1p, 1m, 2p, 2m) to obtain the new selection response.
        Returns: Source dictionary masked with the unsheared mask and with the mean response updated.
        """
        photoz_masks = [(source['bpz_mean'][i] > zlim_low) & (source['bpz_mean'][i] < zlim_high) for i in range(5)]
        source_bin = {}
        source_bin['ra'] = source['ra'][photoz_masks[0]]
        source_bin['dec'] = source['dec'][photoz_masks[0]]
        source_bin['e1'] = source['e1'][photoz_masks[0]]
        source_bin['e2'] = source['e2'][photoz_masks[0]]
        source_bin['psf_e1'] = source['psf_e1'][photoz_masks[0]]
        source_bin['psf_e2'] = source['psf_e2'][photoz_masks[0]]
        source_bin['snr'] = source['snr'][photoz_masks[0]]
        source_bin['size'] = source['size'][photoz_masks[0]]
        source_bin['bpz_mean'] = source['bpz_mean'][0][photoz_masks[0]]
        source_bin['bpz_zmc'] = source['bpz_zmc'][0][photoz_masks[0]]

        if 'v1' in self.config['mastercat_v']:
            R11, _, _ = calibrator.calibrate('e1', mask=photoz_masks)
            R22, _, _ = calibrator.calibrate('e2', mask=photoz_masks)
        if 'v2' in self.config['mastercat_v']:
            R11, _, _ = calibrator.calibrate('e_1', mask=photoz_masks)
            R22, _, _ = calibrator.calibrate('e_2', mask=photoz_masks)
        source_bin['Rmean'] = np.mean([R11, R22])

        return source_bin


    def get_lens(self, lens):
        """
        Given a lens sample, returns ra, dec, jk and weight, in case it exists.
        """

        ra_l = lens['ra']
        dec_l = lens['dec']
        jk_l = lens['jk']
        try:
            w_l = lens['w']
            print 'Weights found in lens catalog.'
        except:
            print 'There are no identified weights for the lenses.'
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
            print 'Weights found in source catalog.'
        except:
            print 'There are no identified weights for the sources.'
            w = np.ones(len(ra_s))

        return ra_s, dec_s, w

    def run_treecorr_jackknife(self, lens, source, type_corr):
        """
        Function that runs treecorr for a given lens and source sample,
        and a given configuration, and paralalizes the measurements
        in different jackknife patches for the lenses and randoms.
        Returns the measurements for each jackknife patch.
        type_corr: string, type of correlation, i.e. NG, NN.
        NG for gammat, NN for wtheta.
        """
        assert type_corr == 'NG' or type_corr == 'NN', 'This type_corr of correlation is not accepted by this function.'

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

        def run_jki(jk):
            """
            Function we use for mutiprocessing.
            jk: Region we run in each core.
            """

            ra_l_jk = ra_l[jk_l == jk]
            dec_l_jk = dec_l[jk_l == jk]
            w_l_jk = w_l[jk_l == jk]

            if type_corr == 'NG':
                if jk == 0: print 'Doing NG correlation.'
                corr = treecorr.NGCorrelation(nbins=self.config['nthbins'], min_sep=self.config['thlims'][0],
                                              max_sep=self.config['thlims'][1], sep_units='arcmin',
                                              bin_slop=self.config['bslop'])

            if type_corr == 'NN':
                if jk == 0: print 'Doing NN correlation.'
                corr = treecorr.NNCorrelation(nbins=self.config['nthbins'], min_sep=self.config['thlims'][0],
                                              max_sep=self.config['thlims'][1], sep_units='arcmin',
                                              bin_slop=self.config['bslop'])

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

                print 'jk, nlens, nsource = ', jk, len(ra_l_jk), len(ra_s[bool_s])
                cat_l = treecorr.Catalog(ra=ra_l_jk, dec=dec_l_jk, w=w_l_jk, ra_units='deg', dec_units='deg')
                if mode == 'data' or mode=='data_y1sources':
                    cat_s = treecorr.Catalog(ra=ra_s[bool_s], dec=dec_s[bool_s], g1=e1[bool_s], g2=e2[bool_s], w=w[bool_s],
                                             ra_units='deg', dec_units='deg')
                if mode == 'mice':
                    cat_s = treecorr.Catalog(ra=ra_s[bool_s], dec=dec_s[bool_s], g1=-e1[bool_s], g2=e2[bool_s], w=w[bool_s],
                                             ra_units='deg', dec_units='deg')
                corr.process(cat_l, cat_s)

                if jk == 0: theta.append(np.exp(corr.logr))
                if type_corr == 'NG':
                    gts[jk].append(corr.xi)
                    gxs[jk].append(corr.xi_im)
                    errs[jk].append(np.sqrt(np.abs(corr.varxi)))
                weights[jk].append(corr.weight)
                npairs[jk].append(corr.npairs)

            else:
                if jk == 0: theta.append(np.exp(corr.logr))
                zeros = np.zeros(self.config['nthbins'])
                if type_corr == 'NG':
                    gts[jk].append(zeros)
                    gxs[jk].append(zeros)
                    errs[jk].append(zeros)
                weights[jk].append(zeros)
                npairs[jk].append(zeros)

        ra_l, dec_l, jk_l, w_l = self.get_lens(lens)
        ra_s, dec_s, w = self.get_source(source)
        e1 = source['e1']
        e2 = source['e2']

        nside = 8
        theta = (90.0 - dec_s) * np.pi / 180.
        phi = ra_s * np.pi / 180.
        pix = hp.ang2pix(nside, theta, phi)

        manager = Manager()
        theta = manager.list()
        weights = [manager.list() for x in range(self.config['njk'])]
        npairs = [manager.list() for x in range(self.config['njk'])]
        gts = [manager.list() for x in range(self.config['njk'])]
        gxs = [manager.list() for x in range(self.config['njk'])]
        errs = [manager.list() for x in range(self.config['njk'])]

        p = mp.Pool(10, worker_init)
        p.map(run_jki, range(self.config['njk']))
        p.close()

        def reshape_manager(obj):
            return (np.array(list(obj))).reshape(self.config['njk'], self.config['nthbins'])

        print 'Reshaping manager... theta'
        theta = (np.array(list(theta))).reshape(self.config['nthbins'])
        print 'Reshaping manager... weights'
        weights = reshape_manager(weights)
        print 'Reshaping manager... npairs'
        npairs = reshape_manager(npairs)

        if type_corr == 'NG':
            print 'Reshaping manager... gts, gx, err'
            gts = reshape_manager(gts)
            gxs = reshape_manager(gxs)
            errs = reshape_manager(errs)
            return theta, gts, gxs, errs, weights, npairs

        if type_corr == 'NN':
            print 'returning NN'
            return theta, weights, npairs

    def run_nk(self, lens, source):
        """
        Uses TreeCorr to compute the NK correlation between lens and source.
        Used to compute scale dependece responses.
        Returns theta and R_nk.
        """

        ra_l, dec_l, jk_l, w_l = self.get_lens(lens)
        ra_s, dec_s, w = self.get_source(source)
        Rgamma = source['Rgamma']

        nk = treecorr.NKCorrelation(nbins=self.config['nthbins'], min_sep=self.config['thlims'][0],
                                    max_sep=self.config['thlims'][1], sep_units='arcmin', bin_slop=self.config['bslop'])

        cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l, w=w_l, ra_units='deg', dec_units='deg')
        cat_s = treecorr.Catalog(ra=ra_s, dec=dec_s, k=Rgamma, w=w, ra_units='deg', dec_units='deg')
        nk.process(cat_l, cat_s)
        theta = np.exp(nk.logr)
        R_nk = nk.xi

        return theta, R_nk

    def compute_Rs(self, e_ix):
        """
        Computes R_s.
        e_ix: Dictionary per each component 1p, 1m, 2p, 2m.
        It can be averaged over all angular scales, or averaged in angular bins using NK correlation.
        Note: In Y1, e_ix is already divided by dgamma in another script, that's why we don't do it here, but might need to be changed for Y3.
        """
        Rs11_mean = e_ix['1p'] - e_ix['1m']
        Rs22_mean = e_ix['2p'] - e_ix['2m']
        Rs_mean = 0.5 * (Rs11_mean + Rs22_mean)
        return Rs_mean

    def load_responses_mean(self, sbin):
        """
        Loads responses for sbin.
        Note: sbin can be 'notomo', which is the union of all the bins.
        Returns R_mean
        """
        path_save = os.path.join(self.paths['runs'], 'metacal_%s' % self.config['metacal_v'],
                                 'zbinsource_%s' % self.config['zbinsource_v']) + '/'
        R_mean, _, _ = np.loadtxt(path_save + 'responses_mean_%s' % sbin, unpack=True)
        return R_mean

    def save_responses_mean(self, Rgamma_mean, Rs_mean, sbin):
        """
        Puts together Rgamma and Rs and saves them.
        Returns R_total.
        """
        R_mean = Rgamma_mean + Rs_mean
        responses = np.array([R_mean, Rgamma_mean, Rs_mean])
        path_save = os.path.join(self.paths['runs'], 'metacal_%s' % self.config['metacal_v'],
                                 'zbinsource_%s' % self.config['zbinsource_v']) + '/'
        np.savetxt(path_save + 'responses_mean_%s' % sbin, responses.reshape(1, responses.shape[0]),
                   header='R_mean, Rgamma_mean, Rs_mean')
        return R_mean

    def run_responses_nk_tomo(self, lens, source, sbin):
        """
        Function that computes scale dependent responses for each lens-source combination.
        Uses NK TreeCorr correlation.
        """

        print 'Be aware: Y1 mode of responses. Not dividing by dgamma here.'
        R_nk_ix = {}  # Scale dependent R for component i,x for a given selection s, divided by Delta gamma.
        # x: p, m
        components = ['1p', '1m', '2p', '2m']
        for comp in components:
            source_selection = pf.getdata(
                self.paths['y1'] + 'metacal_sel_responses/metacal_sel_responses_sa%s_%s.fits' % (sbin[1], comp))

            source_component_ix = {
                'ra': source_selection['ra'],
                'dec': source_selection['dec'],
                'Rgamma': source_selection['Riisx']
            }

            theta, R_nk_ix[comp] = self.run_nk(lens, source_component_ix)

        theta, Rgamma_nk = self.run_nk(lens, source)
        Rs_nk = self.compute_Rs(R_nk_ix)

        R_nk = Rgamma_nk + Rs_nk
        return theta, R_nk, Rgamma_nk, Rs_nk

    def run_responses_mean_notomo(self, Rgamma):
        """
        Rgamma: (R11 + R22)/2 for each galaxy.
        Averages Rgamma in combination of all source bins and computes and averages R_s too.
        Then it computes the mean of R_total.
        Saves R_total, Rgamma, R_s in file and returns R_total.
        """

        print 'Be aware: Y1 mode of responses. Not dividing by dgamma here.'
        e_ix = {}  # ellipticities for component i for a given selection s, divided by Delta gamma.
        # x: p, m
        components = ['1p', '1m', '2p', '2m']
        for comp in components:

            print comp
            e_ix_allbins = np.zeros(0)
            for sbin in zbins['sbins']:
                # Appending responses source bin sbin.
                source_selection = pf.getdata(
                    self.paths['y1'] + 'metacal_sel_responses/metacal_sel_responses_sa%s_%s.fits' % (sbin[1], comp))
                e_ix_allbins = np.append(e_ix_allbins, source_selection['Riisx'])

            e_ix[comp] = np.mean(e_ix_allbins)

        Rgamma_mean = np.mean(Rgamma)
        Rs_mean = self.compute_Rs(e_ix)
        R_mean = self.save_responses_mean(Rgamma_mean, Rs_mean, 'notomo')

        return R_mean

    def save_runs(self, path_test, theta, gts, gxs, errs, weights, npairs, random_bool):
        """
        Function to save the measurements in each jk patch.
        Currently not used, just useful for debugging sometimes.
        """
        rp = '_rp' if random_bool else ''
        np.savetxt(path_test + 'theta' + rp, theta, header='theta [arcmin]')
        np.savetxt(path_test + 'gts' + rp, gts, header='gts for all jackknife regions')
        np.savetxt(path_test + 'gxs' + rp, gxs, header='gxs for all jackknife regions')
        np.savetxt(path_test + 'errs' + rp, errs, header='errs for all jackknife regions')
        np.savetxt(path_test + 'weights' + rp, weights, header='weights for all jackknife regions')
        np.savetxt(path_test + 'npairs' + rp, npairs, header='npairs for all jackknife regions')

    def compute_boost_factor(self, jk_l, jk_r, wnum, wnum_r):
        """
        Computes the boost factor for a given set of weights for the lenses and randoms.
        Uses Eq. from Sheldon et al. (2004)
        It does so for each N-1 jk regions.
        """

        bf_jk = []
        jks = np.arange(self.config['njk'])
        ratio_n_jk = np.zeros(self.config['njk'])
        for jk in jks:
            ratio_n_jk[jk] = float((jk_r != jk).sum()) / (jk_l != jk).sum()
            bf_jk.append(ratio_n_jk[jk] * wnum[jk] / wnum_r[jk])

        bf_jk = np.array(bf_jk)
        return bf_jk

    def numerators_jackknife(self, gts, gxs, ws):
        """
        Given a path and a filename it loads the corresponding file for th, gt, gx, weights, etc.
        Each file containts a quantity for all jk regions, together in the same file.
        Note: Assumes the run was done with process, not process_cross treecorr functions.
        (Advantage of process vs process_cross is that process allows you to save the shape noise err)
        Returns:
        Theta, and the numerators of gt, gx and weights, to be combined later.
        """
        # Obtain the numerators
        gts = gts * ws
        gxs = gxs * ws

        # Construct array taking the sum of all gt from each jk patch except for one patch, different each time
        gt_num = np.array([np.sum(np.delete(np.array(gts), i, axis=0), axis=0) for i in range(len(gts))])
        gx_num = np.array([np.sum(np.delete(np.array(gxs), i, axis=0), axis=0) for i in range(len(gts))])
        w_num = np.array([np.sum(np.delete(np.array(ws), i, axis=0), axis=0) for i in range(len(gts))])

        return gt_num, gx_num, w_num

    def process_run(self, all, theta, path_test, end):
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


class Measurement(GGL):
    """
    Subclass that runs the gglensing measurement for all the lens-source bin combinations.
    Includes:
    - Mean response calculation.
    - Random points subtraction.
    - Jackknife covariance calculation.
    - Boost factors calculation.
    """

    def __init__(self, config, paths, zbins, plotting):
        GGL.__init__(self, config, paths)
        self.zbins = zbins
        self.plotting = plotting

    def get_path_test(self, lbin, sbin):
        return os.path.join(self.paths['runs_config'], 'measurement', lbin + '_' + sbin) + '/'

    def get_path_test_allzbins(self):
        return os.path.join(self.paths['runs_config'], 'measurement') + '/'

    def get_twopointfile_name(self):
        return os.path.join(self.get_path_test_allzbins() + 'gammat_twopointfile.fits')

    def run(self):
        if mode == 'data':
            lens_all = pf.getdata(self.paths['lens'])
            random_all = pf.getdata(self.paths['randoms'])
            source_all, calibrator = self.load_metacal()

        if mode == 'data_y1sources':
            lens_all = pf.getdata(self.paths['lens'])
            random_all = pf.getdata(self.paths['randoms'])

        if mode == 'mice':
            lens_all = pf.getdata(self.paths['lens_mice'])
            random_all = pf.getdata(self.paths['randoms_mice'])
            source_all = pf.getdata(self.paths['source_mice'])

        for sbin in zbins['sbins']:
    		print 'Running measurement for source %s.' % sbin

		if mode == 'data':
		    source = self.load_metacal_bin(source_all, calibrator, zlim_low=zbins[sbin][0], zlim_high=zbins[sbin][1])
		    R = source['Rmean']

		if mode == 'data_y1sources':
		    source = pf.getdata(self.paths['y1'] + 'metacal_sel_sa%s.fits'%sbin[1])

    		if mode == 'mice':
    		    """
    		    In this case there are no responses, so we set it to one.
    		    """
    		    R = 1.
                    source = source_all[(source_all['z'] > zbins[sbin][0]) & (source_all['z'] < zbins[sbin][1])]

    		for l, lbin in enumerate(zbins['lbins']):
    		    print 'Running measurement for lens %s.' % lbin
    		    path_test = self.get_path_test(lbin, sbin)
    		    make_directory(path_test)

    		    lens = lens_all[(lens_all['z'] > zbins[lbin][0]) & (lens_all['z'] < zbins[lbin][1])]

    		    theta, gts, gxs, errs, weights, npairs = self.run_treecorr_jackknife(lens, source, 'NG')
    		    self.save_runs(path_test, theta, gts, gxs, errs, weights, npairs, False)
    		    gtnum, gxnum, wnum = self.numerators_jackknife(gts, gxs, weights)

    		    if mode == 'data':
    			random = random_all[(random_all['z'] > zbins[lbin][0]) & (random_all['z'] < zbins[lbin][1])]
    		    if mode == 'mice':
    			random = random_all[l*len(random_all)/len(zbins['lbins']):(l+1)*len(random_all)/len(zbins['lbins'])]

    		    theta, gts, gxs, errs, weights, npairs = self.run_treecorr_jackknife(random, source, 'NG')
    		    self.save_runs(path_test, theta, gts, gxs, errs, weights, npairs, True)
    		    gtnum_r, gxnum_r, wnum_r = self.numerators_jackknife(gts, gxs, weights)

    		    gt_all = (gtnum / wnum) / R - (gtnum_r / wnum_r) / R
    		    gx_all = (gxnum / wnum) / R - (gxnum_r / wnum_r) / R

    		    bf_all = self.compute_boost_factor(lens['jk'], random['jk'], wnum, wnum_r)

    		    self.process_run(gt_all, theta, path_test, 'gt')
    		    self.process_run(gx_all, theta, path_test, 'gx')
    		    self.process_run((gtnum_r / wnum_r) / R, theta, path_test, 'randoms')
    		    self.process_run(bf_all, theta, path_test, 'boost_factor')

    def save_gammat_2pointfile(self):
        """
        Save the gammat measurements, N(z)'s and jackknife covariance into the 2point format file.
        Creates a:
        - SpectrumMeasurement obejct: In which the gammat measurements are saved. 
        - Kernel object: The N(z)'s are here. 
        - CovarianceMatrixInfo object: In which the jackknife covariance is saved.

        Then, it builds the TwoPointFile objet from the above objects,
        and saves it to a file.
        """

        gt_length = self.config['nthbins'] * len(self.zbins['lbins']) * len(self.zbins['sbins'])
        gt_values = np.zeros(gt_length, dtype=float)
        bin1 = np.zeros(gt_length, dtype=int)
        bin2 = np.zeros_like(bin1)
        angular_bin = np.zeros_like(bin1)
        angle = np.zeros_like(gt_values)
        dv_start = 0
        cov = np.zeros((gt_length, gt_length))

        for l in range(0, len(zbins['lbins'])):
            for s in range(len(zbins['sbins'])):
                path_test = self.get_path_test(zbins['lbins'][l], zbins['sbins'][s])
                theta, gt, gt_err = np.loadtxt(path_test + 'mean_gt', unpack=True)
                cov_ls = np.loadtxt(path_test + 'cov_gt')

                bin_pair_inds = np.arange(dv_start, dv_start + self.config['nthbins'])
                gt_values[bin_pair_inds] = gt
                bin1[bin_pair_inds] = l+1
                bin2[bin_pair_inds] = s+1
                angular_bin[bin_pair_inds] = np.arange(self.config['nthbins'])
                angle[bin_pair_inds] = theta
                dv_start += self.config['nthbins']
                cov[bin_pair_inds[0]:bin_pair_inds[-1]+1, bin_pair_inds] = cov_ls

        # Load Y1 twopoint file to get the N(z)'s for the blinding script
        y1 = twopoint.TwoPointFile.from_fits('y1_2pt_NG_mcal_1110.fits')
        y1_lensnz = y1.get_kernel('nz_lens')
        y1_sourcenz = y1.get_kernel('nz_source')

        gammat = twopoint.SpectrumMeasurement('gammat', (bin1, bin2),
                                                     (twopoint.Types.galaxy_position_real,
                                                      twopoint.Types.galaxy_shear_plus_real),
                                                     ('nz_lens', 'nz_source'), 'SAMPLE', angular_bin, gt_values,
                                                     angle=angle, angle_unit='arcmin')

        cov_mat_info = twopoint.CovarianceMatrixInfo('COVMAT', ['gammat'], [gt_length], cov)

        gammat_twopoint = twopoint.TwoPointFile([gammat], [y1_lensnz, y1_sourcenz], windows=None, covmat_info=cov_mat_info)

        twopointfile_unblind = self.get_twopointfile_name()

        # Remove file if it exists already because to_fits function doesn't overwrite
        if os.path.isfile(twopointfile_unblind):
            os.system('rm %s'%(twopointfile_unblind))

        gammat_twopoint.to_fits(twopointfile_unblind)



    def save_boostfactors_2pointfile(self):
        """
        Save the boost factors into the 2point format file.
        """

        bf_length = self.config['nthbins'] * len(self.zbins['lbins']) * len(self.zbins['sbins'])
        bf_values = np.zeros(bf_length, dtype=float)
        bin1 = np.zeros(bf_length, dtype=int)
        bin2 = np.zeros_like(bin1)
        angular_bin = np.zeros_like(bin1)
        angle = np.zeros_like(bf_values)
        dv_start = 0

        for l in range(0, len(zbins['lbins'])):
            for s in range(len(zbins['sbins'])):
                path_test = self.get_path_test(zbins['lbins'][l], zbins['sbins'][s])
                theta, bf, bf_err = np.loadtxt(path_test + 'mean_boost_factor', unpack=True)

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

    def plot(self):
        """"
        Makes plot of the fiducial measurement for all redshift bins.
        """

        cmap = self.plotting['cmap']
        cmap_step = 0.25
        title_source = self.plotting['catname']

        # Figure
        fig, ax = plt.subplots(2, 3, figsize=(10, 6), sharey=False, sharex=False, gridspec_kw={'height_ratios': [1, 1]})
        fig.subplots_adjust(hspace=0.0, wspace=0.00)

        for l in range(0, len(zbins['lbins'])):

            # To iterate between the three columns and two lines
            j = 0 if l < 3 else 1
            ax[j][l % 3].axvspan(2.5, self.plotting['th_limit'][l], color='gray', alpha=0.2)

            for s in range(len(zbins['sbins'])):

                path_test = self.get_path_test(zbins['lbins'][l], zbins['sbins'][s])
                if os.path.isfile(path_test + 'mean_gt'):
                    th, gt, err = np.loadtxt(path_test + 'mean_gt', unpack=True)

                    mask_neg = gt < 0
                    mask_pos = gt > 0

                    chi2, ndf = self.get_chi2(path_test, 'gt')
                    ax[j][l % 3].errorbar(th[mask_neg] * (1 + 0.05 * s), -gt[mask_neg], err[mask_neg], fmt='.', mfc='None',
                                          mec=plt.get_cmap(cmap)(cmap_step * s), ecolor=plt.get_cmap(cmap)(cmap_step * s), capsize=2)
                    ax[j][l % 3].errorbar(th[mask_pos] * (1 + 0.05 * s), gt[mask_pos], err[mask_pos], fmt='.',
                                          color=plt.get_cmap(cmap)(cmap_step * s),
                                          mec=plt.get_cmap(cmap)(cmap_step * s), label=self.plotting['redshift_s'][s], capsize=2)

                    ax[j][l % 3].set_xlim(2.5, 300)
                    ax[j][l % 3].set_ylim(10 ** (-6), 10 ** (-2))
                    ax[j][l % 3].set_xscale('log')
                    ax[j][l % 3].set_yscale('log')

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

        ax[1][0].set_ylim(10 ** (-6), 0.999 * 10 ** (-2))
        ax[1][1].set_ylim(10 ** (-6), 0.999 * 10 ** (-2))
        # handles, labels = ax[0][0].get_legend_handles_labels()
        fig.delaxes(ax[1, 2])
        # ax[1][1].legend(handles[::-1], labels[::-1], frameon=True, fancybox = True,prop={'size':12}, numpoints = 1, loc='center left', bbox_to_anchor=(1, 0.5))
        ax[0][0].legend(frameon=False, fancybox=True, prop={'size': 12}, numpoints=1, loc='center',
                        bbox_to_anchor=(2.45, -0.52))
        fig.suptitle(title_source, fontsize=16)
        fig.subplots_adjust(top=0.93)
        self.save_plot('plot_measurement')



    def plot_from_twopointfile(self):
        
        """"
        Makes plot of the fiducial measurement for all redshift bins, from a twopoint file, like Y1 style.
        It also uses the twopoint plotting functions to plot the measurements in a different style,
        the covariance and the N(z)'s.
        Useful to plot the blinded measurements (now the default). 
        """

        filename = self.get_twopointfile_name()
        gammat_file = twopoint.TwoPointFile.from_fits('%s_BLINDED.fits'%filename[:-5])

        gammat = gammat_file.spectra[0]
        pairs = gammat.bin_pairs
        npairs = len(pairs)
        # It starts with 1, not with 0
        bins_l = np.transpose(pairs)[0]
        bins_s = np.transpose(pairs)[1]
        nbins_l = np.max(bins_l)
        nbins_s = np.max(bins_s)
        assert len(zbins['lbins']) == nbins_l, 'Number of lens bins in info does not match with the one in the two-point file.'
        assert len(zbins['sbins']) == nbins_s, 'Number of source bins in info does not match with the one in the two-point file.'

        cmap = self.plotting['cmap']
        cmap_step = 0.25
        title_source = self.plotting['catname']

        # Figure
        fig, ax = plt.subplots(2, 3, figsize=(10, 6), sharey=False, sharex=False, gridspec_kw={'height_ratios': [1, 1]})
        fig.subplots_adjust(hspace=0.0, wspace=0.00)

        for l in range(0, len(zbins['lbins'])):

            # To iterate between the three columns and two lines
            j = 0 if l < 3 else 1
            ax[j][l % 3].axvspan(2.5, self.plotting['th_limit'][l], color='gray', alpha=0.2)

            for s in range(len(zbins['sbins'])):

                    path_test = self.get_path_test(zbins['lbins'][l], zbins['sbins'][s])
                    th, gt = gammat.get_pair(l+1, s+1)
                    err = gammat.get_error(l+1, s+1)

                    mask_neg = gt < 0
                    mask_pos = gt > 0

                    chi2, ndf = self.get_chi2(path_test, 'gt')
                    ax[j][l % 3].errorbar(th[mask_neg] * (1 + 0.05 * s), -gt[mask_neg], err[mask_neg], fmt='.', mfc='None',
                                          mec=plt.get_cmap(cmap)(cmap_step * s), ecolor=plt.get_cmap(cmap)(cmap_step * s), capsize=2)
                    ax[j][l % 3].errorbar(th[mask_pos] * (1 + 0.05 * s), gt[mask_pos], err[mask_pos], fmt='.',
                                          color=plt.get_cmap(cmap)(cmap_step * s),
                                          mec=plt.get_cmap(cmap)(cmap_step * s), label=self.plotting['redshift_s'][s], capsize=2)

                    ax[j][l % 3].set_xlim(2.5, 300)
                    ax[j][l % 3].set_ylim(10 ** (-6), 10 ** (-2))
                    ax[j][l % 3].set_xscale('log')
                    ax[j][l % 3].set_yscale('log')

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

        ax[1][0].set_ylim(10 ** (-6), 0.999 * 10 ** (-2))
        ax[1][1].set_ylim(10 ** (-6), 0.999 * 10 ** (-2))
        # handles, labels = ax[0][0].get_legend_handles_labels()
        fig.delaxes(ax[1, 2])
        # ax[1][1].legend(handles[::-1], labels[::-1], frameon=True, fancybox = True,prop={'size':12}, numpoints = 1, loc='center left', bbox_to_anchor=(1, 0.5))
        ax[0][0].legend(frameon=False, fancybox=True, prop={'size': 12}, numpoints=1, loc='center',
                        bbox_to_anchor=(2.45, -0.52))
        fig.suptitle(title_source, fontsize=16)
        fig.subplots_adjust(top=0.93)
        self.save_plot('plot_measurement_BLINDED')

        # Use twopoint library to make the rest of the plots
        gammat_file.plots(self.paths['plots_config'] + 'gammat_twopointfile_BLINDED', blind_yaxis=True, latex = False)


    def plot_boostfactors(self):

        cmap = self.plotting['cmap']
        fig, ax = plt.subplots(4, 5, figsize=(12.5, 9.375), sharey=True, sharex=True)
        fig.subplots_adjust(hspace=0.1, wspace=0.1)
        c1 = plt.get_cmap(cmap)(0.)

        for l in range(0, len(zbins['lbins'])):
            for s in range(len(zbins['sbins'])):

                path_test = self.get_path_test(zbins['lbins'][l], zbins['sbins'][s])
                theta, bf, bf_err = np.loadtxt(path_test + 'mean_boost_factor', unpack=True)

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
                    ax[s][l].set_ylabel('%s\n' % plotting['redshift_s'][s] + r'Boost factors', size='larger',
                                        linespacing=3)
                if s == 0:
                    ax[s][l].set_title(plotting['redshift_l'][l], size='larger')

                ax[s][l].axvspan(2.5, self.plotting['th_limit'][l], color='gray', alpha=0.2)

        ax[0][4].legend(frameon=False, fontsize=16, loc='lower right')
        self.save_plot('boost_factors')

    def plot_randoms(self):
        """
        Makes plot of the tangential shear around random points.
        """

        labels = [plotting['catname']]
        c = 0  # If adding im3shape, s=1
        markers = ['o', '^']
        fs = 18  # fontsize
        cmap = plotting['cmap']
        cmap_step = 0.25
        fig, ax = plt.subplots(4, 5, figsize=(17.25, 13.8), sharey=True, sharex=True)
        fig.subplots_adjust(hspace=0.0, wspace=0.0)

        for l in range(0, len(zbins['lbins'])):

            for s in range(len(zbins['sbins'])):

                path_test = self.get_path_test(zbins['lbins'][l], zbins['sbins'][s])
                th, gt, err = np.loadtxt(path_test + 'mean_randoms', unpack=True)

                ax[s][l].axhline(y=0, ls=':', color='k')
                ax[s][l].errorbar(th * (1 + 0.07 * c), gt, err, fmt=markers[c], color=plt.get_cmap(cmap)(cmap_step * s),
                                  mec=plt.get_cmap(cmap)(cmap_step * s), label=labels[c], capsize=1.3)

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

        ax[0][4].legend(frameon=False, fancybox=True, prop={'size': 13}, numpoints=1, loc='upper right')
        self.save_plot('plot_randoms')

    def plot_gammax(self):
        """
        Makes plot of the cross-component gammax.
        Top panel: Plot of the gammax measurement for a single lens-source combination.
        Bottom panel: chi2 distribution from all lens-source combinations.
        """

        labels = [plotting['catname']]
        c = 0  # for metacal
        markers = ['o', '^']
        fs = 12  # fontsize
        cmap = plotting['cmap']
        cmap_step = 0.25
        c1 = plt.get_cmap(cmap)(cmap_step * 1.5)
        c2 = plt.get_cmap(cmap)(cmap_step * 3)
        colors = [c1, c2]
        fig, ax = plt.subplots(2, 1, figsize=(5, 8), sharey=False, sharex=False)
        fig.subplots_adjust(hspace=0.2, wspace=0.)

        # TOP panel
        l = 0
        s = 0
        path_test = self.get_path_test(zbins['lbins'][l], zbins['sbins'][s])
        th, gx, err = np.loadtxt(path_test + 'mean_gx', unpack=True)
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
        for l in range(0, len(zbins['lbins'])):
            for s in range(len(zbins['sbins'])):
                path_test = self.get_path_test(zbins['lbins'][l], zbins['sbins'][s])
                chi2, ndf = self.get_chi2(path_test, 'gx')
                chi2s.append(chi2)
        chi2s = np.array(chi2s)
        print len(chi2s)
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


class Responses(GGL):
    """
    Subclass that obtains the scale dependence responses (NK correlations) for all the lens-source bin combinations. Both for randoms and lenses.
    """

    def __init__(self, config, paths, zbins, plotting):
        GGL.__init__(self, config, paths)
        self.zbins = zbins
        self.plotting = plotting

    def get_path_test(self, lbin, sbin):
        return os.path.join(self.paths['runs_config'], 'responses_nk', lbin + '_' + sbin) + '/'

    def save_responses_nk(self, path_test, responses, end):
        np.savetxt(path_test + 'responses_nk_%s' % end, responses, header='theta(arcmin), R_nk, Rgamma_nk, Rs_nk')

    def load_responses_nk(self, path_test, end):
        theta, R_nk, Rgamma_nk, Rs_nk = np.loadtxt(path_test + 'responses_nk_%s' % end, unpack=True)
        return theta, R_nk, Rgamma_nk, Rs_nk

    def run_tomo_nk(self):
        """
        Runs the NK responses between lenses and sources.
        Runs for lenses and randoms too.
        """

        lens_all = pf.getdata(self.paths['y1'] + 'lens.fits')
        random_all = pf.getdata(self.paths['y1'] + 'random.fits')

        for sbin in zbins['sbins']:

            print 'Running measurement for source %s.' % sbin
            source = pf.getdata(self.paths['y1'] + 'metacal_sel_sa%s.fits' % sbin[1])

            for lbin in zbins['lbins']:
                print 'Running measurement for lens %s.' % lbin
                path_test = self.get_path_test(lbin, sbin)
                make_directory(path_test)

                lens = lens_all[(lens_all['z'] > zbins[lbin][0]) & (lens_all['z'] < zbins[lbin][1])]
                theta, R_nk, Rgamma_nk, Rs_nk = self.run_responses_nk_tomo(lens, source, sbin)
                self.save_responses_nk(path_test, zip(theta, R_nk, Rgamma_nk, Rs_nk), 'lens')

                random = random_all[(random_all['z'] > zbins[lbin][0]) & (random_all['z'] < zbins[lbin][1])]
                theta, R_nk, Rgamma_nk, Rs_nk = self.run_responses_nk_tomo(random, source, sbin)
                self.save_responses_nk(path_test, zip(theta, R_nk, Rgamma_nk, Rs_nk), 'random')

    def plot(self, lens_random):
        """
        Makes plot comparing the NK responses to the mean ones.
        lens_random: string, can be lens or random.
        Indicates which is the foreground sample when computing the NK correlations.
        """

        cmap = plotting['cmap']
        cmap_step = 0.25
        c1 = plt.get_cmap(cmap)(0.)
        c2 = plt.get_cmap(cmap)(0.6)
        fig, ax = plt.subplots(4, 5, figsize=(16.5, 13.2), sharey='row', sharex=True)
        fig.subplots_adjust(hspace=0.1, wspace=0.1)

        for l in range(0, len(zbins['lbins'])):

            for s in range(len(zbins['sbins'])):

                path_test = self.get_path_test(zbins['lbins'][l], zbins['sbins'][s])
                theta, R_nk, _, _ = self.load_responses_nk(path_test, lens_random)
                R_mean = self.load_responses_mean(zbins['sbins'][s])

                ax[s][l].plot(theta, [R_mean] * len(theta), '-', lw=2, color=c2, mec=c2, label=r'$R_{\mathrm{mean}}$')
                ax[s][l].plot(theta, R_nk, '-', lw=2, color=c1, mec=c1, label=r'$R_{\mathrm{nk}}$')
                ax[s][l].plot(theta, [np.mean(R_nk)] * len(theta), ':', lw=2, color=c1, mec=c1,
                              label=r'$\overline{R_{\mathrm{nk}}}$')
                ax[s][l].set_xscale('log')
                ax[s][l].set_xlim(self.config['thlims'][0], self.config['thlims'][1])
                ax[0][l].set_ylim(0.724, 0.729)
                ax[2][l].set_ylim(0.633, 0.638)
                ax[s][l].xaxis.set_major_formatter(ticker.FormatStrFormatter('$%0.0f$'))
                ax[s][l].tick_params(axis='both', which='major', labelsize='larger')
                ax[s][l].tick_params(axis='both', which='minor', labelsize='larger')

                if s == 3:
                    ax[s][l].set_xlabel(r'$\theta$ [arcmin]', size='larger')
                if l == 0:
                    ax[s][l].set_ylabel('%s\n' % plotting['redshift_s'][s] + r'Responses', size='larger', linespacing=3)
                if s == 0:
                    ax[s][l].set_title(plotting['redshift_l'][l], size='larger')

                ax[s][l].text(0.5, 0.85,
                              r'$\Delta R/R = %0.2f \%%$' % (100 * np.mean((R_nk - R_mean) / (R_mean + R_nk) * 2)),
                              horizontalalignment='center', verticalalignment='center', transform=ax[s][l].transAxes,
                              fontsize='larger')

        ax[0][4].legend(frameon=False, fontsize=16, loc='lower right')
        self.save_plot('plot_responses_%s' % lens_random)


class TestStars(GGL):
    """
    SubClass to test if the tangential shear around stars is consistent with zero.
    Uses no tomography for the source sample.
    """

    def __init__(self, config, paths, zbins, plotting):
        GGL.__init__(self, config, paths)
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


class TestPSF(GGL):
    """
    SubClass to test if the psf residuals are compatible with zero.
    Uses no tomography for the lens sample.
    """

    def __init__(self, config, paths, zbins, plotting):
        GGL.__init__(self, config, paths)
        self.zbins = zbins
        self.plotting = plotting

    def get_path_test(self):
        return os.path.join(self.paths['runs_config'], 'psfresiduals') + '/'

    def save_psf_residuals(self, ra_lims, dec_lims):
        """
        Computes the psf residuals for the r band and saves them to a file.
        """

        info = pf.open(paths['y1base'] + 'cats/y1a1-v13/exposure_info_y1a1-v13.fits')[1].data
        exp = info['exp']
        print 'Original info:', len(exp)
        fil = info['filter']
        flag_info = info['flag']

        min_ra, max_ra = ra_lims
        min_dec, max_dec = dec_lims

        # Use only r-band
        mask_info = (fil == 'r')
        exp = exp[mask_info]
        print 'Only r band:', len(exp)

        # To not have repeated exposure names due to different ccds (in info file every ccd has a line)
        exp_unique = np.unique(exp)
        print 'All ccds:', len(exp_unique)
        print exp_unique

        ra_all, dec_all, psf1_all, psf2_all, res1_all, res2_all, mag_all = [], [], [], [], [], [], []

        for i in range(len(exp_unique)):
            # Load every exposure
            data = pf.open(paths['y1base'] + 'cats/y1a1-v13/psf_cats/%s_exppsf.fits' % exp_unique[i])[1].data
            ra = data['ra']
            dec = data['dec']
            e1 = data['e1']
            e2 = data['e2']
            psf1 = data['psf_e1']
            psf2 = data['psf_e2']
            # Resiual psf is the difference between the measurement of the psf (e1) and the model of the psf(psfex)
            # at the position of the stars, that is the only place you can measure the psf
            res1 = e1 - psf1
            res2 = e2 - psf2
            flag = data['flag']
            mag = data['mag']
            # Use only reserved stars and spt
            mask = ((dec > min_dec) & (dec < max_dec) & (flag == 64))
            ra_all.extend(ra[mask])
            dec_all.extend(dec[mask])
            psf1_all.extend(psf1[mask])
            psf2_all.extend(psf2[mask])
            res1_all.extend(res1[mask])
            res2_all.extend(res2[mask])
            mag_all.extend(mag[mask])
        # Convert to numpy array and simplify name
        # plt.hist(mag_all, bins=25)
        ra = np.array(ra_all)
        dec = np.array(dec_all)
        psf1 = np.array(psf1_all)
        psf2 = np.array(psf2_all)
        res1 = np.array(res1_all)
        res2 = np.array(res2_all)
        # w = np.ones(len(ra))
        print 'Number of exposures:', len(ra)

        c1 = pf.Column(name='RA', format='E', array=ra)
        c2 = pf.Column(name='DEC', format='E', array=dec)
        c3 = pf.Column(name='E1', format='E', array=res1)
        c4 = pf.Column(name='E2', format='E', array=res2)

        CC = [c1, c2, c3, c4]
        hdu = pf.new_table(CC, nrows=len(ra))
        hdu.writeto(self.paths['y1'] + 'psfresiduals.fits', clobber=True)

    def save_psf_residuals_y3(self, ra_lims, dec_lims):
        """
        Computes the psf residuals for the r band and saves them to a file.
        """

        # info = pf.open(paths['y1base'] + 'cats/y1a1-v13/exposure_info_y1a1-v13.fits')[1].data
        # exp = info['exp']
        # print 'Original info:', len(exp)
        # fil = info['filter']
        # flag_info = info['flag']
        #
        min_ra, max_ra = ra_lims
        min_dec, max_dec = dec_lims
        #
        # # Use only r-band
        # mask_info = (fil == 'r')
        # exp = exp[mask_info]
        # print 'Only r band:', len(exp)
        #
        # # To not have repeated exposure names due to different ccds (in info file every ccd has a line)
        # exp_unique = np.unique(exp)
        # print 'All ccds:', len(exp_unique)
        # print exp_unique

        ra_all_bandr, dec_all_bandr, psf1_all_bandr, psf2_all_bandr, res1_all_bandr, res2_all_bandr, mag_all_bandr, exp_all_bandr, band_all_bandr = [], [], [], [], [], [], [], [], []

        ra_all_bandi, dec_all_bandi, psf1_all_bandi, psf2_all_bandi, res1_all_bandi, res2_all_bandi, mag_all_bandi, exp_all_bandi, band_all_bandi = [], [], [], [], [], [], [], [], []

        ra_all_bandz, dec_all_bandz, psf1_all_bandz, psf2_all_bandz, res1_all_bandz, res2_all_bandz, mag_all_bandz, exp_all_bandz, band_all_bandz = [], [], [], [], [], [], [], [], []

        exposures = []

        y3_exp_dir = self.paths['y3_exp']

        for root, directories, filenames in os.walk(y3_exp_dir):
            for directory in directories:
                exposures.append(directory)

        exposures_int = np.array([int(exp) for exp in exposures])
        exp_unique = np.unique(exposures_int)

        print len(exp_unique)

        for i in range(len(exp_unique)):

            if np.mod(i,100) == 0:
                print i

            # Load every exposure
            data = pf.open(y3_exp_dir + str(exp_unique[i]) + '/exp_psf_cat_%s.fits' % exp_unique[i])[1].data
            info = pf.open(y3_exp_dir + str(exp_unique[i]) + '/exp_psf_cat_%s.fits' % exp_unique[i])[2].data

            flag = np.array(info['flag'])

            if np.sum(flag) == 0:
                band = info['band'][0][0]
                ra = data['ra']
                dec = data['dec']
                obs_e1 = data['obs_e1']
                obs_e2 = data['obs_e2']
                piff_e1 = data['piff_e1']
                piff_e2 = data['piff_e2']
                # Resiual psf is the difference between the measurement of the psf (e1) and the model of the psf(psfex)
                # at the position of the stars, that is the only place you can measure the psf
                res1 = obs_e1 - piff_e1
                res2 = obs_e2 - piff_e2
                obs_flag = data['obs_flag']
                mag = data['mag']
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

                elif band == 'i':
                    ra_all_bandi.extend(ra[mask])
                    dec_all_bandi.extend(dec[mask])
                    psf1_all_bandi.extend(piff_e1[mask])
                    psf2_all_bandi.extend(piff_e2[mask])
                    res1_all_bandi.extend(res1[mask])
                    res2_all_bandi.extend(res2[mask])
                    mag_all_bandi.extend(mag[mask])

                elif band == 'z':
                    ra_all_bandz.extend(ra[mask])
                    dec_all_bandz.extend(dec[mask])
                    psf1_all_bandz.extend(piff_e1[mask])
                    psf2_all_bandz.extend(piff_e2[mask])
                    res1_all_bandz.extend(res1[mask])
                    res2_all_bandz.extend(res2[mask])
                    mag_all_bandz.extend(mag[mask])

                else:
                    print 'no correct band alloted'

        # Convert to numpy array and simplify name
        # plt.hist(mag_all, bins=25)
        ra_bandr = np.array(ra_all_bandr)
        dec_bandr = np.array(dec_all_bandr)
        psf1_bandr = np.array(psf1_all_bandr)
        psf2_bandr = np.array(psf2_all_bandr)
        res1_bandr = np.array(res1_all_bandr)
        res2_bandr = np.array(res2_all_bandr)

        print 'Number of exposures:', len(ra_bandr)

        c1 = pf.Column(name='RA', format='E', array=ra_bandr)
        c2 = pf.Column(name='DEC', format='E', array=dec_bandr)
        c3 = pf.Column(name='E1', format='E', array=res1_bandr)
        c4 = pf.Column(name='E2', format='E', array=res2_bandr)

        CC = [c1, c2, c3, c4]
        hdu = pf.new_table(CC, nrows=len(ra_bandr))
        hdu.writeto(self.paths['y3'] + 'psfresiduals_bandr.fits', clobber=True)

        ra_bandi = np.array(ra_all_bandi)
        dec_bandi = np.array(dec_all_bandi)
        psf1_bandi = np.array(psf1_all_bandi)
        psf2_bandi = np.array(psf2_all_bandi)
        res1_bandi = np.array(res1_all_bandi)
        res2_bandi = np.array(res2_all_bandi)

        print 'Number of exposures bandi:', len(ra_bandi)

        c1 = pf.Column(name='RA', format='E', array=ra_bandi)
        c2 = pf.Column(name='DEC', format='E', array=dec_bandi)
        c3 = pf.Column(name='E1', format='E', array=res1_bandi)
        c4 = pf.Column(name='E2', format='E', array=res2_bandi)

        CC = [c1, c2, c3, c4]
        hdu = pf.new_table(CC, nrows=len(ra_bandi))
        hdu.writeto(self.paths['y3'] + 'psfresiduals_bandi.fits', clobber=True)

        ra_bandz = np.array(ra_all_bandz)
        dec_bandz = np.array(dec_all_bandz)
        psf1_bandz = np.array(psf1_all_bandz)
        psf2_bandz = np.array(psf2_all_bandz)
        res1_bandz = np.array(res1_all_bandz)
        res2_bandz = np.array(res2_all_bandz)

        # w = np.ones(len(ra))
        print 'Number of exposures bandz:', len(ra_bandz)

        c1 = pf.Column(name='RA', format='E', array=ra_bandz)
        c2 = pf.Column(name='DEC', format='E', array=dec_bandz)
        c3 = pf.Column(name='E1', format='E', array=res1_bandz)
        c4 = pf.Column(name='E2', format='E', array=res2_bandz)

        CC = [c1, c2, c3, c4]
        hdu = pf.new_table(CC, nrows=len(ra_bandz))
        hdu.writeto(self.paths['y3'] + 'psfresiduals_bandz.fits', clobber=True)

    def run(self):
        """
        Obtains the tangential component of the psf residuals around lenses, with random point subtraction.
        Obtains the corresponding jackknife covariance.
        """

        lens = pf.getdata(self.paths['y1'] + 'lens.fits')
        masklens = ((lens['z'] > zbins['l1'][0]) & (lens['z'] < zbins['l5'][1]))
        lens = lens[masklens]

        random = pf.getdata(self.paths['y1'] + 'random.fits')
        maskrandom = ((random['z'] > zbins['l1'][0]) & (random['z'] < zbins['l5'][1]))
        print zbins['l1'][0]
        print zbins['l5'][1]
        random = random[maskrandom]

        psfres = pf.getdata(self.paths['y1'] + 'psfresiduals.fits')

        path_test = self.get_path_test()
        make_directory(path_test)

        print 'PSF residuals around lenses...'
        theta, gts, gxs, errs, weights, npairs = self.run_treecorr_jackknife(lens, psfres, 'NG')
        gtnum, gxnum, wnum = self.numerators_jackknife(gts, gxs, weights)

        print 'PSF residuals around randoms...'
        theta, gts, gxs, errs, weights, npairs = self.run_treecorr_jackknife(random, psfres, 'NG')
        gtnum_r, gxnum_r, wnum_r = self.numerators_jackknife(gts, gxs, weights)

        gt_all = gtnum / wnum - gtnum_r / wnum_r

        self.process_run(gt_all, theta, path_test, 'gt')

    def run_y3(self,bands):
        """
        Obtains the tangential component of the psf residuals around lenses, with random point subtraction.
        Obtains the corresponding jackknife covariance.
        """

        lens = pf.getdata(self.paths['y3'] + 'lens.fits')
        masklens = ((lens['z'] > zbins['l1'][0]) & (lens['z'] < zbins['l5'][1]))
        lens = lens[masklens]

        random = pf.getdata(self.paths['y3'] + 'random.fits')
        maskrandom = ((random['z'] > zbins['l1'][0]) & (random['z'] < zbins['l5'][1]))
        print zbins['l1'][0]
        print zbins['l5'][1]
        random = random[maskrandom]

        for band in bands:

            psfres = pf.getdata(self.paths['y3'] + 'psfresiduals_band' + band + '.fits')

            path_test = self.get_path_test()
            make_directory(path_test)

            print 'PSF residuals around lenses...'
            theta, gts, gxs, errs, weights, npairs = self.run_treecorr_jackknife(lens, psfres, 'NG')
            gtnum, gxnum, wnum = self.numerators_jackknife(gts, gxs, weights)

            print 'PSF residuals around randoms...'
            theta, gts, gxs, errs, weights, npairs = self.run_treecorr_jackknife(random, psfres, 'NG')
            gtnum_r, gxnum_r, wnum_r = self.numerators_jackknife(gts, gxs, weights)

            gt_all = gtnum / wnum - gtnum_r / wnum_r

            self.process_run(gt_all, theta, path_test, 'gt')

    def plot(self):
        """
        Makes plot of the psf resdiuals.
        """
        cmap = plotting['cmap']
        c1 = plt.get_cmap(cmap)(0)
        titles_l = r'$0.15 < z < 0.90 $'
        title_redmagic = 'redMaGiC'
        fig, ax = plt.subplots(1, 1, figsize=(4, 4))
        fig.subplots_adjust(hspace=0.0, wspace=0.00)

        path_test = self.get_path_test()
        th, gt, err = np.loadtxt(path_test + 'mean_gt', unpack=True)
        ax.errorbar(th, gt, err, fmt='.', color=c1, mec=c1, markersize=5.7, capsize=1.4)
        ax.set_xlim(2.5, 250)
        ax.set_ylim(-1.8 * 10 ** (-5), 1.8 * 10 ** (-5))
        ax.set_xscale('log')
        ax.axhline(y=0, ls=':', color='k')
        ax.text(0.5, 0.85, titles_l, horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(0.5, 0.92, title_redmagic, horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, fontsize=12)

        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%d$'))

        ax.set_xlabel(r'$\theta$ [arcmin]', size='large')
        ax.set_ylabel(r'$\gamma_{t,\mathrm{PSF\, residuals}}$', size='large')

        chi2, ndf = self.get_chi2(path_test, 'gt')

        ax.text(0.7, 0.18, r'Null $\chi^2$/ndf ',
                horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=10)
        ax.text(0.7, 0.11, r'$%0.1f/%d$' % (chi2, ndf),
                horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=12, color=c1)

        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        self.save_plot('plot_psfresiduals')


class TestSizeSNR(GGL):
    """
    SubClass to test if the tangential shear signal has no dependence on source size or S/N. Using the first lens bin and all sources, as in Y1.
    The variable size_snr on run will be 'size' or 'snr' and will specify the test.
    """

    def __init__(self, config, paths, zbins, plotting, source_nofz_pars):
        GGL.__init__(self, config, paths)
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

    def run_responses_mean_notomo_size_snr(self, Rgamma, size_snr, cut, high_low):
        """
        Computes responses when there is an extra selection on size or snr.
        For all the sources (no tomography).
        - Rgamma: Rgamma for the high or low selection in size or snr.
        - size_snr: string, either size or snr.
        - cut: median of size or snr for the whole source sample.
        - high_low: string, either high or low.
        """

        e_ix = {}  # ellipticities for component i for a given selection s, divided by Delta gamma.
        # x: p, m
        components = ['1p', '1m', '2p', '2m']
        for comp in components:
            print comp
            e_ix_allbins = np.zeros(0)
            par_ix_allbins = np.zeros(0)  # size or snr in the selection ix, for instance, size_1p, size_2p, etc.
            for sbin in zbins['sbins']:
                # Appending responses source bin sbin.
                source_selection = pf.getdata(
                    self.paths['y1'] + 'metacal_sel_responses/metacal_sel_responses_sa%s_%s.fits' % (sbin[1], comp))
                print 'dgamma is now dividing in the function that loads metacal and writes the file on it!'
                e_ix_allbins = np.append(e_ix_allbins, source_selection['Riisx'])
                par_ix_allbins = np.append(par_ix_allbins, source_selection['%s_ix' % size_snr])

            # For each component, mask the high or low part of e_ix
            if high_low == 'high':
                mask = np.where(par_ix_allbins > cut)
            if high_low == 'low':
                mask = np.where(par_ix_allbins <= cut)

            e_ix[comp] = np.mean(e_ix_allbins[mask])

        Rgamma_mean = np.mean(Rgamma)
        Rs_mean = self.compute_Rs(e_ix)
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
        lens_all = pf.getdata(self.paths['y1'] + 'lens.fits')
        lens = lens_all[(lens_all['z'] > zbins['l1'][0]) & (lens_all['z'] < zbins['l1'][1])]
        random_all = pf.getdata(self.paths['y1'] + 'random.fits')
        random = random_all[(random_all['z'] > zbins['l1'][0]) & (random_all['z'] < zbins['l1'][1])]

        sources = pf.getdata(self.paths['y1'] + 'metacal_sel_allbins.fits')

        # Source size or snr split, selecting the two halves of the split
        par = sources[size_snr]
        cut = np.median(par)
        print 'len(sources)', len(sources)
        maskl = par <= cut
        maskh = par > cut

        path_test = self.get_path_test(size_snr)
        make_directory(path_test)
        print size_snr
        print 'NEW len(sources[maskl])', len(sources[maskl])
        print 'NEW len(sources[maskh])', len(sources[maskh])

        # Computing the measurements for the split halves, both around lenses and randoms
        theta, gtsl, gxsl, errsl, weightsl, npairsl = self.run_treecorr_jackknife(lens, sources[maskl], 'NG')
        self.save_runs(path_test, theta, gtsl, gxsl, errsl, weightsl, npairsl, False)
        gtlnum, gxlnum, wlnum = self.numerators_jackknife(gtsl, gxsl, weightsl)

        theta, gtsl_r, gxsl_r, errsl_r, weightsl_r, npairsl_r = self.run_treecorr_jackknife(random, sources[maskl],
                                                                                            'NG')
        self.save_runs(path_test, theta, gtsl_r, gxsl_r, errsl_r, weightsl_r, npairsl_r, True)
        gtlnum_r, gxlnum_r, wlnum_r = self.numerators_jackknife(gtsl_r, gxsl_r, weightsl_r)

        theta, gtsh, gxsh, errsh, weightsh, npairsh = self.run_treecorr_jackknife(lens, sources[maskh], 'NG')
        gthnum, gxhnum, whnum = self.numerators_jackknife(gtsh, gxsh, weightsh)

        theta, gtsh_r, gxsh_r, errsh_r, weightsh_r, npairsh_r = self.run_treecorr_jackknife(random, sources[maskh],
                                                                                            'NG')
        gthnum_r, gxhnum_r, whnum_r = self.numerators_jackknife(gtsh_r, gxsh_r, weightsh_r)

        # Computing the responses for the split halves
        Rl = self.run_responses_mean_notomo_size_snr(sources['Rgamma'][maskl], size_snr, cut, 'low')
        Rh = self.run_responses_mean_notomo_size_snr(sources['Rgamma'][maskh], size_snr, cut, 'high')

        # Combining measurements and responses to get gammat
        gtl_all = (gtlnum / wlnum) / Rl - (gtlnum_r / wlnum_r) / Rl
        gth_all = (gthnum / whnum) / Rh - (gthnum_r / whnum_r) / Rh
        np.savetxt(path_test + 'gtl_all', gtl_all)
        np.savetxt(path_test + 'gth_all', gth_all)

        # Getting the data ratio using the simulations
        sims, cov_sims = self.load_sims()
        ratio, err_ratio = self.ratio_from_sims(theta, gtl_all, gth_all, sims, cov_sims)

        # Load N(z)'s and corrects the mean using Cosmos calibration.
        zl, nzl, zs, nzsl, nzsh = self.load_nzs(size_snr)

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

    def plot(self):
        cmap = plotting['cmap']
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


class TestSysMaps(GGL):
    """
    SubClass to test if the tangential shear signal has no dependence on observational conditions such as seeing, airmass etc, in each band griz.
    Using the first lens bin and all sources, as in Y1.
    The variable map can be: 'airmass', 'fwhm', 'maglimit', 'skybrite'. We iterate over them.
    The variable band can be: 'g', 'r', 'i', 'z'. We iterate over them. In Y1 we only used r band in the end, because it was used by im3shape.
    """

    def __init__(self, config, paths, zbins, plotting, source_nofz_pars, sysmaps):
        GGL.__init__(self, config, paths)
        self.zbins = zbins
        self.plotting = plotting
        self.source_nofz_pars = source_nofz_pars
        self.sysmaps = sysmaps

    def get_path_test(self, map, band):
        return os.path.join(self.paths['runs_config'], 'systematics_maps', map, band) + '/'

    def save_systematics_maps_ratios(self, path_test, result_data, result_theory):
        np.savetxt(path_test + 'data', result_data, header='ratio err_ratio')
        np.savetxt(path_test + 'theory', result_theory, header='ratio err_ratio')

    def load_systematics_maps_ratios(self, path_test):
        data, data_err = np.loadtxt(path_test + 'data', unpack=True)
        theory, theory_err = np.loadtxt(path_test + 'theory', unpack=True)
        return data, data_err, theory, theory_err

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

    def load_systematics_map(self, map, band, nside, nested_bool):
        '''
        Loads the systematics map, splits into high and low parts, and plots all maps.
        nested_bool: True if nest, False if ring.
        Returns: pixels corresponding the low half and the high half, for each map.
        '''
        path = os.path.join(self.paths['y1base'], 'cats', 'systematics_maps') + '/'
        sys_map = pf.getdata(
            path + 'Y1A1NEW_COADD_SPT_band_%s/' % band + 'Y1A1NEW_COADD_SPT_band_%s_nside%s_oversamp4_' % (
                band, nside) + self.sysmaps[map] + '.fits')
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

    def load_nzs(self, map, band):
        """
	Loads redshift distributions for lenses, sources and source splits, low(l) and high(h).
	"""
        zl, nzl = np.loadtxt(self.paths['nz_lens'] + 'lens', unpack=True, usecols=(0, 1))
        zs, nzsl, nzsh = np.loadtxt(
            self.paths['y1base'] + 'runs/test_mcal_bpzmof_unblind/nofzs/source_%s_%s_notomo' % (map, band), unpack=True)
        return zl, nzl, zs, nzsl, nzsh

    def radec_to_thetaphi(self, ra, dec):
        """
        Converts ra and dec in degrees to theta and phi.
        Returns theta and phi in radians.
        """
        theta = (90. - dec) * np.pi / 180.
        phi = ra * np.pi / 180.
        return theta, phi

    def run(self, maps, bands):
        """
        Runs gglensing measurment for all maps and bands.
        """

        lens_all = pf.getdata(self.paths['y1'] + 'lens.fits')
        lens = lens_all[(lens_all['z'] > zbins['l1'][0]) & (lens_all['z'] < zbins['l1'][1])]
        random_all = pf.getdata(self.paths['y1'] + 'random.fits')
        random = random_all[(random_all['z'] > zbins['l1'][0]) & (random_all['z'] < zbins['l1'][1])]
        sources = pf.getdata(self.paths['y1'] + 'metacal_sel_allbins.fits')
        R = self.run_responses_mean_notomo(sources['Rgamma'])

        for map in maps:
            print 'Running map %s...' % map
            for band in bands:
                print 'Band %s' % band

                path_test = self.get_path_test(map, band)
                make_directory(path_test)

                # Load and split the systematics map
                pix_low, pix_hi = self.load_systematics_map(map, band, self.sysmaps['nside'],
                                                            self.sysmaps['nested_bool'])

                # Building lenses masks, low and high
                theta_l, phi_l = self.radec_to_thetaphi(lens['ra'], lens['dec'])
                pix_all_l = hp.pixelfunc.ang2pix(self.sysmaps['nside'], theta_l, phi_l, self.sysmaps['nested_bool'])
                maskl_low = np.in1d(pix_all_l, pix_low)
                maskl_hi = np.in1d(pix_all_l, pix_hi)

                # Building randoms masks, low and high
                theta_r, phi_r = self.radec_to_thetaphi(random['ra'], random['dec'])
                pix_all_r = hp.pixelfunc.ang2pix(self.sysmaps['nside'], theta_r, phi_r, self.sysmaps['nested_bool'])
                maskr_low = np.in1d(pix_all_r, pix_low)
                maskr_hi = np.in1d(pix_all_r, pix_hi)

                # Building sources masks, low and high
                theta_s, phi_s = self.radec_to_thetaphi(sources['ra'], sources['dec'])
                pix_all_s = hp.pixelfunc.ang2pix(self.sysmaps['nside'], theta_s, phi_s, self.sysmaps['nested_bool'])
                masks_low = np.in1d(pix_all_s, pix_low)
                masks_hi = np.in1d(pix_all_s, pix_hi)

                # Computing the measurements for the split halves, both around lenses and randoms
                print 'Lenses, low.'
                theta, gtsl, gxsl, errsl, weightsl, npairsl = self.run_treecorr_jackknife(lens[maskl_low],
                                                                                          sources[masks_low], 'NG')
                gtlnum, gxlnum, wlnum = self.numerators_jackknife(gtsl, gxsl, weightsl)

                print 'Randoms, low.'
                theta, gtsl_r, gxsl_r, errsl_r, weightsl_r, npairsl_r = self.run_treecorr_jackknife(random[maskr_low],
                                                                                                    sources[masks_low],
                                                                                                    'NG')
                gtlnum_r, gxlnum_r, wlnum_r = self.numerators_jackknife(gtsl_r, gxsl_r, weightsl_r)

                print 'Lenses, high.'
                theta, gtsh, gxsh, errsh, weightsh, npairsh = self.run_treecorr_jackknife(lens[maskl_hi],
                                                                                          sources[masks_hi], 'NG')
                gthnum, gxhnum, whnum = self.numerators_jackknife(gtsh, gxsh, weightsh)

                print 'Randoms, high.'
                theta, gtsh_r, gxsh_r, errsh_r, weightsh_r, npairsh_r = self.run_treecorr_jackknife(random[maskr_hi],
                                                                                                    sources[masks_hi],
                                                                                                    'NG')
                gthnum_r, gxhnum_r, whnum_r = self.numerators_jackknife(gtsh_r, gxsh_r, weightsh_r)

                # Combining measurements and responses to get gammat
                gtl_all = (gtlnum / wlnum) / R - (gtlnum_r / wlnum_r) / R
                gth_all = (gthnum / whnum) / R - (gthnum_r / whnum_r) / R

                # Getting the data ratio using the simulations
                sims, cov_sims = self.load_sims()
                ratio, err_ratio = self.ratio_from_sims(theta, gtl_all, gth_all, sims, cov_sims)

                # Load N(z)'s
                zl, nzl, zs, nzsl, nzsh = self.load_nzs(map, band)

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
                self.save_systematics_maps_ratios(path_test, result_data, result_theory)

    def plot(self):
        """
        Plots measurments compared to expected theory.
        Set up for four maps, and a single band.
        """
        #labels_cshort = r'\textsc{Metacal}'
        labels_cshort = 'Metacal'
        fontsize = 16
        maps = ['airmass', 'fwhm', 'maglimit', 'skybrite']
        band = 'r'
        cmap = self.plotting['cmap']
        cmap_step = 0.25
        c1 = plt.get_cmap(cmap)(cmap_step * 1.5)
        c2 = plt.get_cmap(cmap)(cmap_step * 3)

        fig, ax = plt.subplots(2, 2, figsize=(6, 6), sharey=True, sharex=True)
        fig.subplots_adjust(hspace=0.0, wspace=0.00)

        for m, map in enumerate(maps):

            path_test = self.get_path_test(map, band)
            data, data_err, theory, theory_err = self.load_systematics_maps_ratios(path_test)
            if m % 2 == 0: c = 0
            if m % 2 == 1: c = 1
            if m == 0 or m == 1: f = 0
            if m == 2 or m == 3: f = 1

            color, x = c1, 1

            ax[f][c].errorbar(x, data, data_err, fmt='o', color=color, mec=color,
                              label='%s' % (self.plotting['catname']))
            x_lin = [x - 0.2, x, x + 0.2]
            ax[f][c].fill_between(x_lin, np.array([theory - theory_err for i in x_lin]),
                                  np.array([theory + theory_err for i in x_lin]), alpha=0.4, edgecolor=color,
                                  facecolor=color,
                                  label=r'$\Sigma_{\mathrm{crit,eff}}^{-1,\,\mathrm{high}}/ \Sigma_{\mathrm{crit,eff}}^{-1,\,\mathrm{low}} %s$' % (
                                      labels_cshort))

            ax[f][c].set_xlim(0.2, 2.8)
            ax[f][c].set_xticklabels([])
            ax[f][c].set_xticks([])
            ax[f][c].set_ylim(0.5, 1.5)
            ax[f][0].set_ylabel(r'$\gamma^{\mathrm{high}}_t /\ \gamma^{\mathrm{low}}_t $', fontsize=fontsize)
            ax[f][c].text(0.5, 0.9, r'\textsc{%s}' % map,
                          horizontalalignment='center', verticalalignment='center', transform=ax[f][c].transAxes,
                          fontsize=12, color='k')

        ax[1][0].legend(frameon=False, loc='best', numpoints=1, fontsize=9.8)
        self.save_plot('systematics_maps')


T = True
F = False

run_measurement = T
run_responses_nk = F
run_stars = F
run_psf = F
run_size_snr = F
run_sysmaps = F


if run_measurement:
    print 'Starting measurement class...'
    gglensing = GGL(config, paths)
    measurement = Measurement(config, paths, zbins, plotting)
    if not plot_blinded:
        measurement.run()
        #measurement.save_boostfactors_2pointfile() #there is a bug here now
        measurement.save_gammat_2pointfile()
        if not blind:
            measurement.plot()
        measurement.plot_boostfactors()
        measurement.plot_randoms()
        measurement.plot_gammax()

    if blind and plot_blinded:
        measurement.plot_from_twopointfile()

if run_responses_nk:
    responses = Responses(config, paths, zbins, plotting)
    responses.run_tomo_nk()
    responses.plot('lens')
    responses.plot('random')

if run_stars:
    stars = TestStars(config, paths, zbins, plotting)
    stars.run('bright')
    stars.run('faint')
    stars.plot()

if run_psf:
    psf = TestPSF(config, paths, zbins, plotting)
    ra_lims = (-1, 361)
    dec_lims = (-90, -35)
    psf.save_psf_residuals_y3(ra_lims, dec_lims)
    # psf.run()
    # psf.plot()

if run_size_snr:
    size_snr = TestSizeSNR(config, paths, zbins, plotting, source_nofz_pars)
    size_snr.run('size')
    size_snr.run('snr')
    size_snr.plot()

if run_sysmaps:
    sysmaps = TestSysMaps(config, paths, zbins, plotting, source_nofz_pars, sysmaps)
    sysmaps.run(['airmass', 'fwhm', 'maglimit', 'skybrite'], ['r'])
    sysmaps.plot()
