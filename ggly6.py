"""
-----------------------------
Functions to run GGL for Y6KP
-----------------------------
"""

import os
import sys
import numpy as np
import gc
from mpi4py import MPI
from astropy.io import fits

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


class GGL(object):

    def __init__(self, input_dir=None):

        " Load the parameters/settings file "
        if input_dir is not None:
            if os.path.exists(input_dir):
                sys.path.append(input_dir)
            else:
                errmsg = '!!!Error: The input directory {} does not exist'.format(input_dir)
                raise Exception(errmsg)
        else:
            errmsg = '!!!Error: Please provide the path to the input directory'
            raise Exception(errmsg)

        import params as par

        self.par = par

        " Import the setup class "
        import setup

        self.ggl_setup = setup.GGL_setup(input_dir=input_dir)

        return
    

    def setup_run(self, data_response=None,
                  lens_bin=None, source_bin=None,
                  zl_lims=None, zs_lims=None):
        """
        Setup the parameters to run code by reading files
        """

        # load lens data
        print('Reading lens data for redshift bin {:d} [{:.2f},{:.2f}] from {}'.format(
                lens_bin+1, zl_lims[0], zl_lims[1], self.par.data_lens))
        (self.ra_l, self.dec_l, self.weight_lens) = self.ggl_setup.load_lens_Y6_maglim(self.par.data_lens, lens_bin)
        # self.ra_l, self.dec_l, self.weight_lens = self.ggl_setup.load_lens_Y3_maglim(self.par.data_lens, zl_lims)

        if not self.par.use_LSSweight:
            print('Warning: Will not load LSS weights for lenses, setting them to 1')
            self.weight_lens = np.ones(len(self.ra_l))

        # load source data
        if self.par.source_cat == 'metadetect':

            print('Reading source data for source bin {:d} [{:.2f},{:.2f}] from {}'.format(
                    source_bin+1, zs_lims[0], zs_lims[1], self.par.data_source))
            
            (self.ra_s, self.dec_s,
             self.e1_s, self.e2_s,
             self.R_g, self.w_g) = self.ggl_setup.load_source_metadetect(self.par.data_source, data_response, source_bin)

            if not self.par.use_response:
                print('Warning: Will set metadetect response to 1')
                self.R_g = 1.0

            if not self.par.use_shearweight:
                print('Warning: Will not load metadetect weights, setting them to 1')
                self.w_g = np.ones(len(self.ra_s))

        elif self.par.source_cat == 'bfd':

            print('Reading source data for source bin {:d} [{:.2f},{:.2f}] from {}'.format(
                    source_bin+1, zs_lims[0], zs_lims[1], self.par.data_source_bfd))
            
            (self.ra_s, self.dec_s,
             self.e1_s, self.e2_s,
             self.P, self.Q0, self.Q1, 
             self.R00, self.R01, self.R11,
             self.w_g) = self.ggl_setup.load_source_bfd(self.par.data_source_bfd,
                                                        self.par.data_source_bfd_binning, self.par.data_source_bfd_mask, source_bin)

            if not self.par.use_shearweight:
                print('Warning: Will not load BFD weights, setting them to 1')
                self.w_g = np.ones(len(self.ra_s))

        else:
            print('Specify which source catalog you want to use in the params.py file! Exiting the program')
            sys.exit(0)

        # load random points data
        if (self.par.use_randoms or self.par.use_boosts or self.par.compute_covariance):
            print('Reading random-point data from {}'.format(self.par.data_randoms))
            (self.ra_rand, self.dec_rand) = self.ggl_setup.load_randoms_Y6(self.par.data_randoms, lens_bin)
            #(self.ra_rand, self.dec_rand) = self.ggl_setup.load_randoms_Y3_maglim(self.par.data_randoms, zl_lims)
        else:
            print('Will not load randoms points data, as it is not needed in current run')
            self.ra_rand = None
            self.dec_rand = None
            
        print('Done reading data')

        return
    

    def create_output_dirs(self):

        # output directory for gamma_t
        print('Saving gamma_t measurements in {}'.format(self.par.path_out_gt))
        if not os.path.exists(self.par.path_out_gt):
            os.makedirs(self.par.path_out_gt)

        # output directory for gamma_x
        if not os.path.exists(self.par.path_out_gx):
            os.makedirs(self.par.path_out_gx)

        # output directory for boost factors
        if not os.path.exists(self.par.path_out_boost):
            os.makedirs(self.par.path_out_boost)

        # setup output path for extra info
        if not os.path.exists(self.par.path_out_extra_gt):
            os.makedirs(self.par.path_out_extra_gt)

        if not os.path.exists(self.par.path_out_gt_rand):
            os.makedirs(self.par.path_out_gt_rand)

        if not os.path.exists(self.par.path_out_gx_rand):
            os.makedirs(self.par.path_out_gx_rand)

        # setup output path for gamma_t shot noise variance
        if not os.path.exists(self.par.path_out_shot_gt):
            os.makedirs(self.par.path_out_shot_gt)

        if self.par.compute_covariance == True:

            # setup output path for gamma_t Jackknife covariance
            if not os.path.exists(self.par.path_JK_cov_gt):
                os.makedirs(self.par.path_JK_cov_gt)

            # setup output path for random-point gamma_t Jackknife covariance
            if not os.path.exists(self.par.path_JK_cov_gt_rand):
                os.makedirs(self.par.path_JK_cov_gt_rand)

            # setup output path for gamma_t Jackknife covariance
            if not os.path.exists(self.par.path_JK_cov_gx):
                os.makedirs(self.par.path_JK_cov_gx)

            # setup output path for gamma_t Jackknife covariance
            if not os.path.exists(self.par.path_JK_cov_bf):
                os.makedirs(self.par.path_JK_cov_bf)
                
        if self.par.path_out_gt[-1] == '/':
            self.path_out_gt_final = self.par.path_out_gt[:-1]
        else:
            self.path_out_gt_final = self.par.path_out_gt
            
        if self.par.path_out_gx[-1] == '/':
            self.path_out_gx_final = self.par.path_out_gx[:-1]
        else:
            self.path_out_gx_final = self.par.path_out_gx
            
        if self.par.use_boosts:
            self.path_out_gt_final += '_bf'
            self.path_out_gx_final += '_bf'
        if self.par.use_randoms:
            self.path_out_gt_final += '_rp'
            self.path_out_gx_final += '_rp'
            
        if not os.path.exists(self.path_out_gt_final):
            os.makedirs(self.path_out_gt_final)
                
        if not os.path.exists(self.path_out_gx_final):
            os.makedirs(self.path_out_gx_final)
        return
    

    def save_results(self,
                     lzind, szind,
                     theta_res, 
                     gammat_total, gammat_res, gammat_rand,
                     gammax_total, gammax_res, gammax_rand,
                     shot_noise_gammat,
                     xi_im, xi_im_rand, xi_npairs, xi_npairs_rand, xi_weight, xi_weight_rand,
                     Rg, sum_w_l, sum_w_r, boosts):

        # gamma_t output directory
        gammat_out = self.par.path_out_gt+'/gammat_l{0}_s{1}.txt'.format(lzind+1, szind+1)
        gammax_out = self.par.path_out_gx+'/gammax_l{0}_s{1}.txt'.format(lzind+1, szind+1)
        extra_out = self.par.path_out_extra_gt+'/gammat_extra_l{0}_s{1}.txt'.format(lzind+1, szind+1)
        randoms_gt_out = self.par.path_out_gt_rand+'/gammat_rand_l{0}_s{1}.txt'.format(lzind+1, szind+1)
        randoms_gx_out = self.par.path_out_gx_rand+'/gammax_rand_l{0}_s{1}.txt'.format(lzind+1, szind+1)
        extra_rand_out = self.par.path_out_extra_gt+'/gammat_rand_extra_l{0}_s{1}.txt'.format(lzind+1, szind+1)
        boosts_out = self.par.path_out_boost+'/boost_l{0}_s{1}.txt'.format(lzind+1, szind+1)
        shot_gammat_out = self.par.path_out_shot_gt+'/shot_noise_gammat_l{0}_s{1}.txt'.format(lzind+1, szind+1)

        #  save results in file
        np.savetxt(gammat_out, np.c_[theta_res, gammat_res], header='theta, gamma_t')
        np.savetxt(gammax_out, np.c_[theta_res, gammax_res], header='theta, gamma_x')
        np.savetxt(extra_out, np.c_[xi_im, xi_npairs, xi_weight, Rg * np.ones(len(theta_res))],
                   header='xi_im, xi_npair, xi_weight, Rg')
        
        np.savetxt(randoms_gt_out, np.c_[theta_res, gammat_rand], header='theta, gamma_t')
        np.savetxt(randoms_gx_out, np.c_[theta_res, xi_im_rand / Rg], header='theta, gamma_x')
        np.savetxt(extra_rand_out, np.c_[xi_im_rand, xi_npairs_rand, xi_weight_rand, 
                                         Rg * np.ones(len(theta_res)), 
                                         sum_w_l * np.ones(len(theta_res)), sum_w_r * np.ones(len(theta_res))],
                   header='xi_im, xi_npair, xi_weight, Rg, sum_w_l, sum_w_r')

        np.savetxt(boosts_out, np.c_[theta_res, boosts], header='theta, boost')
        np.savetxt(shot_gammat_out, np.c_[theta_res, shot_noise_gammat], header='theta, gammat_shot_noise')

        if (self.par.use_boosts or self.par.use_randoms):
            gammat_total_out = self.path_out_gt_final+'/gammat_l{0}_s{1}.txt'.format(lzind+1, szind+1)
            gammax_total_out = self.path_out_gx_final+'/gammax_l{0}_s{1}.txt'.format(lzind+1, szind+1)
            
            np.savetxt(gammat_total_out, np.c_[theta_res, gammat_total], header='theta, gamma_t')
            np.savetxt(gammax_total_out, np.c_[theta_res, gammax_total], header='theta, gamma_t')
        else:
            if not np.all(gammat_total / gammat_res == 1.0):
                errmsg = '!!!Something is wrong, no boost or randoms, but final gammat is not equal to the basic gammat measurement'
                raise Exception(errmsg)
            if not np.all(gammax_total / (xi_im / Rg) == 1.0):
                errmsg = '!!!Something is wrong, no randoms, but final gammax is not equal to the basic gammax measurement'
                raise Exception(errmsg)

        print('Results saved in {}'.format(self.par.out_main_gt))

        return

    
    def save_cov(self, lzind, szind,
                 cov_gammat, cov_gammax, cov_boost, cov_gammat_rand):
        
        cov_gammat_out = (self.par.path_JK_cov_gt+'/cov_gammat_l{0}_s{1}.txt'.format(lzind+1, szind+1))
        cov_gammat_rand_out = (self.par.path_JK_cov_gt_rand+'/cov_gammat_rand_l{0}_s{1}.txt'.format(lzind+1, szind+1))
        cov_gammax_out = (self.par.path_JK_cov_gx+'/cov_gammax_l{0}_s{1}.txt'.format(lzind+1, szind+1))
        cov_boosts_out = (self.par.path_JK_cov_bf+'/cov_boost_l{0}_s{1}.txt'.format(lzind+1, szind+1))
        err_gammat_out = (self.par.path_JK_cov_gt+'/err_gammat_l{0}_s{1}.txt'.format(lzind+1, szind+1))
        
        with open(cov_gammat_out,'wb') as f:
            for line in cov_gammat.T:
                np.savetxt(f, [line])
                
        gt_err = np.sqrt(np.diag(cov_gammat))
        np.savetxt(err_gammat_out, np.c_[theta_res,gt_err], header='theta, gammat_err')
        
        with open(cov_boosts_out,'wb') as f:
            for line in cov_boost.T:
                np.savetxt(f, [line])
                
        with open(cov_gammax_out,'wb') as f:
            for line in cov_gammax.T:
                np.savetxt(f, [line])
                
        with open(cov_gammat_rand_out,'wb') as f:
            for line in cov_gammat_rand.T:
                np.savetxt(f, [line])
                        
        print('Covariances saved in {}'.format(self.path_out_JK))

        return
    

    def save_results_2pt_file(self):
        
        theta_all = []
        gammat_all = []
        
        for l in self.par.l_bins:
            for s in self.par.s_bins:
                gamma_t = np.loadtxt(self.path_out_gt_final+'/gammat_l{0}_s{1}.txt'.format(l+1, s+1))
                theta_all.append(gamma_t[:, 0])
                gammat_all.append(gamma_t[:, 1])
        
        theta_all = np.concatenate(theta_all)                
        gammat_all = np.concatenate(gammat_all)

        dv = fits.open(self.par.dv_input)
        dv[4].data['ANG'] = theta_all
        dv[4].data['VALUE'] = gammat_all
        dv.writeto(self.par.dv_output, overwrite=True)
        
        print('2pt file saved in {}'.format(self.par.dv_output))
        
        return
    

    def save_params(self):
        
        params_dict = {key: value for key, value in vars(self.par).items()
            if not key.startswith('__') and not callable(value)}
        
        with open(self.par.out_main+'/params_content.txt', 'w') as f:
            for key, value in params_dict.items():
                f.write(f'{key} = {value}\n')
                
        print('Parameters saved in {}'.format(self.par.out_main+'/params_content.txt'))
        
        return
    

    def run_gammat(self, lzind, szind):
        """
        Run code to get gamma_t
        """

        # lens redshift cuts
        zl_min, zl_max = self.par.zl_bins[lzind]
        # source redshift cuts
        zs_min, zs_max = self.par.zs_bins[szind]

        # give feedback on progress
        print('\nLens bin {:d} [{:.2f},{:.2f}] x source bin {:d} [{:.2f},{:.2f}]'.format(
            lzind+1, zl_min, zl_max, szind+1, zs_min, zs_max))
        sys.stdout.flush()

        # load data and setup current bin
        self.setup_run(data_response=self.par.data_response[szind],
                       lens_bin=lzind, source_bin=szind,
                       zl_lims=[zl_min, zl_max], zs_lims=[zs_min, zs_max])

        # random points
        if (self.par.use_randoms or self.par.use_boosts):
            print('Num randoms = {:d}'.format(len(self.ra_rand)))
        else:
            print('Will not use random points')

        if (self.par.source_cat == 'metadetect' or self.par.ggl_bfd_approximate == True):
            
            # parameters to parse to treecorr
            params = [self.e1_s, self.e2_s, self.R_g, self.w_g]

            # get gamma_t for defined parameters
            (theta_res, 
             gammat_total, gammat_res, gammat_rand,
             gammax_total, gammax_res, gammax_rand,
             shot_noise_gammat,
             xi_im, xi_im_rand, xi_npairs, xi_npairs_rand, xi_weight, xi_weight_rand,
             Rg, sum_w_l, sum_w_r, 
             boosts) = self.ggl_setup.get_gammat(self.ra_l, self.dec_l, self.ra_s, self.dec_s, 
                                                                self.ra_rand, self.dec_rand, params=params,
                                                                low_mem=self.par.treecorr_low_mem,
                                                                weights=self.weight_lens,
                                                                use_randoms=self.par.use_randoms,
                                                                use_boosts=self.par.use_boosts)

        else:
            
            # parameters to parse to treecorr
            params = [self.Q0, self.Q1, self.R00, self.R01, self.R11]

            # get gamma_t for defined parameters
            (theta_res, 
             gammat_total, gammat_res, gammat_rand,
             gammax_total, gammax_res, gammax_rand,
             shot_noise_gammat,
             xi_im, xi_im_rand, xi_npairs, xi_npairs_rand, xi_weight, xi_weight_rand,
             Rg, sum_w_l, sum_w_r, 
             boosts) = self.ggl_setup.get_gammat_bfd(self.ra_l, self.dec_l, self.ra_s, self.dec_s,
                                                                    self.ra_rand, self.dec_rand, params=params,
                                                                    low_mem=self.par.treecorr_low_mem,
                                                                    weights=self.weight_lens,
                                                                    use_randoms=self.par.use_randoms,
                                                                    use_boosts=self.par.use_boosts)

        self.save_results(lzind, szind,
                          theta_res, 
                          gammat_total, gammat_res, gammat_rand,
                          gammax_total, gammax_res, gammax_rand,
                          shot_noise_gammat,
                          xi_im, xi_im_rand, xi_npairs, xi_npairs_rand, xi_weight, xi_weight_rand,
                          Rg, sum_w_l, sum_w_r, boosts)

        # clear up memory
        del self.ra_l, self.dec_l, self.ra_s, self.dec_s
        del self.ra_rand, self.dec_rand
        del self.weight_lens, self.w_g
        if (self.par.source_cat == 'metadetect' or self.par.ggl_bfd_approximate == True):
            del self.e1_s, self.e2_s, self.R_g
        else:
            del self.Q0, self.Q1, self.R00, self.R01, self.R11
        
        
        return
    

    def run_gammat_and_cov(self, lzind, szind):
        """
        Run code to get gamma_t and its Jackknife covariance
        """

        # lens redshift cuts
        zl_min, zl_max = self.par.zl_bins[lzind]
        # source redshift cuts
        zs_min, zs_max = self.par.zs_bins[szind]

        # give feedback on progress
        print('\nLens bin {:d} [{:.2f},{:.2f}] x source bin {:d} [{:.2f},{:.2f}]'.format(
            lzind+1, zl_min, zl_max, szind+1, zs_min, zs_max))
        sys.stdout.flush()

        # load data and setup current bin
        self.setup_run(data_response=self.par.data_response[szind],
                       lens_bin=lzind, source_bin=szind,
                       zl_lims=[zl_min, zl_max], zs_lims=[zs_min, zs_max])
        
        # random points
        if (self.par.use_randoms or self.par.use_boosts):
            print('Num randoms = {:d}'.format(len(self.ra_rand)))
        else:
            print('Will not use random points')

        if (self.par.source_cat == 'metadetect' or self.par.ggl_bfd_approximate == True):
            
            # parameters to parse to treecorr
            params = [self.e1_s, self.e2_s, self.R_g, self.w_g]

            # get gamma_t for defined parameters
            (theta_res, 
             gammat_total, gammat_res, gammat_rand,
             gammax_total, gammax_res, gammax_rand,
             cov_gammat, cov_gammax, cov_boost, cov_gammat_rand,
             shot_noise_gammat,
             xi_im, xi_im_rand, xi_npairs, xi_npairs_rand, xi_weight, xi_weight_rand,
             Rg, sum_w_l, sum_w_r, 
             boosts) = self.ggl_setup.get_gammat_and_covariance(self.ra_l, self.dec_l, self.ra_s, self.dec_s, 
                                                                self.ra_rand, self.dec_rand, params=params,
                                                                low_mem=self.par.treecorr_low_mem,
                                                                weights=self.weight_lens,
                                                                use_randoms=self.par.use_randoms,
                                                                use_boosts=self.par.use_boosts)

        else:
            
            # parameters to parse to treecorr
            params = [self.Q0, self.Q1, self.R00, self.R01, self.R11]

            # get gamma_t for defined parameters
            (theta_res, 
             gammat_total, gammat_res, gammat_rand,
             gammax_total, gammax_res, gammax_rand,
             cov_gammat, cov_gammax, cov_boost, cov_gammat_rand,
             shot_noise_gammat,
             xi_im, xi_im_rand, xi_npairs, xi_npairs_rand, xi_weight, xi_weight_rand,
             Rg, sum_w_l, sum_w_r, 
             boosts) = self.ggl_setup.get_gammat_and_covariance_bfd(self.ra_l, self.dec_l, self.ra_s, self.dec_s,
                                                                    self.ra_rand, self.dec_rand, params=params,
                                                                    low_mem=self.par.treecorr_low_mem,
                                                                    weights=self.weight_lens,
                                                                    use_randoms=self.par.use_randoms,
                                                                    use_boosts=self.par.use_boosts)

        self.save_results(lzind, szind,
                          theta_res, 
                          gammat_total, gammat_res, gammat_rand,
                          gammax_total, gammax_res, gammax_rand,
                          shot_noise_gammat,
                          xi_im, xi_im_rand, xi_npairs, xi_npairs_rand, xi_weight, xi_weight_rand,
                          Rg, sum_w_l, sum_w_r, boosts)
        
        self.save_cov(lzind, szind,
                      cov_gammat, cov_gammax, cov_boost, cov_gammat_rand)

        # clear up memory
        del self.ra_l, self.dec_l, self.ra_s, self.dec_s
        del self.ra_rand, self.dec_rand
        del self.weight_lens, self.w_g
        if (self.par.source_cat == 'metadetect' or self.par.ggl_bfd_approximate == True):
            del self.e1_s, self.e2_s, self.R_g
        else:
            del self.Q0, self.Q1, self.R00, self.R01, self.R11

        gc.collect()

        return

    
    def run_gammat_measurement(self):
        """
        Run code to get gamma_t
        
        Options:
        - parallel / not parallel
        - gamma_t only / with its Jackknife covariance using Treecorr
        
        ------
        Results are saved in file
        """

        print('Computing gamma_t with bin slop={:.2f} and resolution={:d}:'.format(
                self.par.bin_slop, self.par.nside))
        print('Running TreeCorr with theta=[{:.1f},{:.1f}] over {:d} angular bins'.format(
                self.par.theta_lims[0], self.par.theta_lims[1], self.par.ang_nbins))
        
        self.create_output_dirs()

        if self.par.run_parallel == True:

            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()

            # Define the total number of lens and source bins
            num_lens_bins = len(self.par.l_bins)
            num_source_bins = len(self.par.s_bins)

            # Calculate the chunks per process
            chunks_per_process = (num_lens_bins * num_source_bins) // size

            # Calculate the start and end chunk for this process
            start_chunk = rank * chunks_per_process
            end_chunk = start_chunk + chunks_per_process

            # Iterate over the assigned chunks
            for chunk in range(start_chunk, end_chunk):
                lzind = chunk // num_source_bins
                szind = chunk % num_source_bins
                
                if self.par.compute_covariance == True:
                    self.run_gammat_and_cov(lzind, szind)
                else:
                    self.run_gammat(lzind, szind)
        else:
            
            for lzind in self.par.l_bins:
                for szind in self.par.s_bins:
                    
                    if self.par.compute_covariance == True:
                        self.run_gammat_and_cov(lzind, szind)
                    else:
                        self.run_gammat(lzind, szind)
                        
        print('\nDone calculating gamma_t')
        
        if((len(self.par.l_bins) == 6) & (len(self.par.s_bins) == 4)):
            self.save_results_2pt_file()

        # Save the content of params.py to a text file
        self.save_params()
           
        print('\nWARNING: all the results saved are UNBLINDED, need to run blinding\n')

        return