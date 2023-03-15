"""
-----------------------------
Functions to run GGL for Y6KP
-----------------------------
"""
import os
import sys
import numpy as np

class GGL(object):

    def __init__(self, input_dir=None):

        " load the parameters/settings file "
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

        " import the setup class "
        import setup
        self.ggl_setup = setup.GGL_setup(input_dir=input_dir)

        return
    
    def setup_run(self, lens_file=None, lens_dir=None, 
                  source_file=None,
                  randoms_file=None, randoms_dir=None, 
                  lens_bin=None, source_bin=None, 
                  zl_lims=None, zs_lims=None,
                  load_sources=True):
        """
        Setup the parameters to run code by reading files
        """
        print( "Setting up things to run code:" )

        # load lens data
        print("Reading lens data for redshift bin %d from %s..."%(lens_bin+1,lens_file))

        # read lens galaxy data
        self.ra_l, self.dec_l, self.z_l = self.self.ggl_setup.load_lens_Y3_maglim(lens_dir+'/'+lens_file)

        # LSS weights
        if self.par.use_LSSweight:
            self.weight_lens = self.self.ggl_setup.get_weightLSS(self.ra_l, self.dec_l, mask_file=self.par.weightLSS_file, 
                                                                 NSIDE=self.par.nside, nest=self.par.weightLSS_nest, key=self.lens_bin_key)
        else:
            self.weight_lens = np.ones(len(self.ra_l))

        # load source data
        if load_sources:
            print("Reading source data for source bin %d from %s..."%(source_bin+1,self.par.data_source))
            # read source galaxy data
            source, source_5sels, source_calibrator = self.self.ggl_setup.load_source_metacal_5sels(source_file)
            # mask sources
            (self.ra_s, self.dec_s, 
             self.e1_s, self.e2_s, 
             self.R_g, self.w_g) = self.self.ggl_setup.mask_source_metacal_5sels(source['ra'], source['dec'], source['e1'], source['e2'], 
                                                                                 source_5sels, source_calibrator, 
                                                                                 zs_bin=source_bin)
            
        # load random points data
        if self.par.use_randoms or self.par.use_boosts:
            self.ra_rand, self.dec_rand = self.self.ggl_setup.load_randoms_Y3_maglim(randoms_dir+'/'+randoms_file)
        else:
            print("No need to load source data, skipping this part...")

        print( "Done reading data" )

        return

    def run_gammat(self):
        """
        Run code to get gamma_t

        output
        ------
        results are saved in file
        """
        # output directory for gamma_t
        path_out_gt = self.par.path_out_gt
        if path_out_gt[-1] is not '/': path_out_gt+='/'
        if not os.path.exists(path_out_gt):
            os.makedirs(path_out_gt)

        # output directory for gamma_x
        path_out_gx = self.par.path_out_gx
        if path_out_gx[-1] is not '/': path_out_gx+='/'
        if not os.path.exists(path_out_gx):
            os.makedirs(path_out_gx)

        # output directory for boost factors
        path_out_boost = self.par.path_out_boost
        if path_out_boost[-1] is not '/': path_out_boost+='/'
        if not os.path.exists(path_out_boost):
            os.makedirs(path_out_boost)

        # setup output path for extra info
        path_out_extra = self.par.path_out_extra_gt
        if path_out_extra[-1] is not '/': path_out_extra+='/'
        if not os.path.exists(path_out_extra):
            os.makedirs(path_out_extra)

        # setup output path for randoms
        path_out_rand = self.par.path_out_rand
        if path_out_rand[-1] is not '/': path_out_rand+='/'
        if not os.path.exists(path_out_rand):
            os.makedirs(path_out_rand)
        #
        path_out_gt_rand = self.par.path_out_gt_rand
        if path_out_gt_rand[-1] is not '/': path_out_gt_rand+='/'
        if not os.path.exists(path_out_gt_rand):
            os.makedirs(path_out_gt_rand)
        #
        path_out_gx_rand = self.par.path_out_gx_rand
        if path_out_gx_rand[-1] is not '/': path_out_gx_rand+='/'
        if not os.path.exists(path_out_gx_rand):
            os.makedirs(path_out_gx_rand)

        # setup output path for gamma_t shot noise variance
        path_out_shot_gt = self.par.path_out_shot_gt
        if path_out_shot_gt[-1] is not '/': path_out_shot_gt+='/'
        if not os.path.exists(path_out_shot_gt):
            os.makedirs(path_out_shot_gt)

        # print feedback
        print( "Working on gamma_t calculation with bin slop=%.3f and resolution=%d:"%(self.par.bin_slop,self.par.nside) )
        print( "Running treecorr with theta=[%.1f,%.1f] over %d angular bins"%(self.par.theta_lims[0],self.par.theta_lims[1],self.par.ang_nbins) )

        # run code to get gamma_t
        for lzind in self.par.l_bins:
            # lens redshift cuts
            zl_min, zl_max = self.par.zl_bins[lzind]

            for szind in self.par.s_bins:
                # source redshift cuts
                zs_min, zs_max = self.par.zs_bins[szind]

                # give feedback on progress
                print( "  Doing: lens bin l%d [%.2f,%.2f] x source bin %d [%.2f,%.2f]"%(lzind+1,zl_min,zl_max,szind+1,zs_min,zs_max) )
                # gamma_t output directory
                gammat_out = path_out_gt+'/gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                gammax_out = path_out_gx+'/gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                extra_out = path_out_extra+'/gammat_extra_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                extra_rand_out = path_out_extra+'/gammat_rand_extra_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                randoms_gt_out = path_out_gt_rand+'/gammat_rand_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                randoms_gx_out = path_out_gx_rand+'/gammax_rand_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                boosts_out = path_out_boost+'/boost_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                shot_gammat_out = path_out_shot_gt+'/shot_noise_gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                
                # load data and setup current bin
                self.setup_run(lens_file=self.par.data_lens[lzind], lens_dir=self.par.lens_dir, 
                               randoms_file=self.par.data_randoms[lzind], randoms_dir=self.par.randoms_dir, 
                               source_file=self.par.data_source[szind], 
                               lens_bin=lzind, source_bin=szind, zl_lims=[zl_min,zl_max], zs_lims=[zs_min,zs_max])

                print('Number of lenses=',len(self.ra_l))

                # random points
                if self.par.use_randoms or self.par.use_boosts:
                    print('Number of randoms=',len(self.ra_rand))
                else:
                    print('Will not use random points')

                # parameters to parse to treecorr
                params = [self.e1_s,self.e2_s,self.R_g,self.w_g]
                
                # get gamma_t for defined parameters
                (theta_res, gammat_total, gammat_res, gammat_rand, 
                 shot_noise_gammat,
                 xi_im, xi_im_rand, xi_npairs, xi_npairs_rand, xi_weight, xi_weight_rand, 
                 Rg, sum_w_l, sum_w_r, 
                 boosts) = self.ggl_setup.get_gammat(self.ra_l, self.dec_l, self.ra_s, self.dec_s,
                                                     ra_rand=self.ra_rand, dec_rand=self.dec_rand, 
                                                     params=params, low_mem=self.par.treecorr_low_mem, weights=self.weight_lens, 
                                                     use_randoms=self.par.use_randoms, use_boosts=self.par.use_boosts)
                # save gamma_x
                np.savetxt(gammax_out, np.c_[theta_res,xi_im/Rg], header='theta, gamma_x')

                # save shot noise
                np.savetxt(shot_gammat_out, np.c_[theta_res,shot_noise_gammat], header='theta, gammat_shot_noise')
                
                # save results in file
                np.savetxt(randoms_gt_out, np.c_[theta_res,gammat_rand], header='theta, gamma_t')
                np.savetxt(randoms_gx_out, np.c_[theta_res,xi_im_rand/Rg], header='theta, gamma_x')
                np.savetxt(extra_rand_out, np.c_[xi_im_rand,xi_npairs_rand,xi_weight_rand,
                                                    Rg*np.ones(len(theta_res)),sum_w_l*np.ones(len(theta_res)),sum_w_r*np.ones(len(theta_res))], header='xi_im, xi_npair, xi_weight, Rg, sum_w_l, sum_w_r')

                # save results in file
                np.savetxt(gammat_out, np.c_[theta_res,gammat_res], header='theta, gamma_t')
                np.savetxt(extra_out, np.c_[xi_im,xi_npairs,xi_weight,Rg*np.ones(len(theta_res))], header='xi_im, xi_npair, xi_weight, Rg')
                np.savetxt(boosts_out, np.c_[theta_res,boosts], header='theta, boost')

                # piece together components to get gamma_t with RP subtraction and/or boost factors applied
                if path_out_gt[-1]=='/':
                    path_out_gt_final = path_out_gt[:-1]
                else:
                    path_out_gt_final = path_out_gt
                if self.par.use_boosts:
                    path_out_gt_final += '_bf'
                if self.par.use_randoms:
                    path_out_gt_final += '_RP'
                path_out_gt_final += '/'
                gammat_total_out = path_out_gt_final+'gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                # setup output path
                if not os.path.exists(path_out_gt_final):
                    os.makedirs(path_out_gt_final)
                np.savetxt(gammat_total_out, np.c_[theta_res,gammat_total], header='theta, gamma_t')

                # piece together components to get gamma_x with RP subtraction and/or boost factors applied
                gammax_total = xi_im/Rg - xi_im_rand/Rg
                if path_out_gx[-1]=='/':
                    path_out_gx_final = path_out_gx[:-1]
                else:
                    path_out_gx_final = path_out_gx
                if self.par.use_randoms:
                    path_out_gx_final += '_RP'
                path_out_gx_final += '/'
                gammax_total_out = path_out_gx_final+'gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                # setup output path
                if not os.path.exists(path_out_gx_final):
                    os.makedirs(path_out_gx_final)
                np.savetxt(gammax_total_out, np.c_[theta_res,gammax_total], header='theta, gamma_x')

                # give feedback on progress
                print( "  Results saved in: %s"%gammat_out )
                print( "--Done\n" )

                # clear up memory to avoid out-of-memory issues
                del self.ra_l, self.dec_l, self.ra_s, self.dec_s
                del ra_rand, dec_rand
                del self.e1_s,self.e2_s,self.R_g,self.w_g
                del self.weight_lens

        print( "Done calculating gamma_t \n" )
        return
    
    def run_gammat_and_cov(self):
        """
        Run code to get gamma_t and its covariance using only Treecorr

        output
        ------
        results are saved in file
        """
        # output directory for gamma_t
        path_out_gt = self.par.path_out_gt
        if path_out_gt[-1] is not '/': path_out_gt+='/'
        if not os.path.exists(path_out_gt):
            os.makedirs(path_out_gt)

        # output directory for gamma_x
        path_out_gx = self.par.path_out_gx
        if path_out_gx[-1] is not '/': path_out_gx+='/'
        if not os.path.exists(path_out_gx):
            os.makedirs(path_out_gx)

        # output directory for boost factors
        path_out_boost = self.par.path_out_boost
        if path_out_boost[-1] is not '/': path_out_boost+='/'
        if not os.path.exists(path_out_boost):
            os.makedirs(path_out_boost)

        # setup output path for extra info
        path_out_extra = self.par.path_out_extra_gt
        if path_out_extra[-1] is not '/': path_out_extra+='/'
        if not os.path.exists(path_out_extra):
            os.makedirs(path_out_extra)

        # setup output path for randoms
        path_out_rand = self.par.path_out_rand
        if path_out_rand[-1] is not '/': path_out_rand+='/'
        if not os.path.exists(path_out_rand):
            os.makedirs(path_out_rand)
        #
        path_out_gt_rand = self.par.path_out_gt_rand
        if path_out_gt_rand[-1] is not '/': path_out_gt_rand+='/'
        if not os.path.exists(path_out_gt_rand):
            os.makedirs(path_out_gt_rand)
        #
        path_out_gx_rand = self.par.path_out_gx_rand
        if path_out_gx_rand[-1] is not '/': path_out_gx_rand+='/'
        if not os.path.exists(path_out_gx_rand):
            os.makedirs(path_out_gx_rand)

        # setup output path for gamma_t Jackknife covariance
        path_JK_cov_gt = self.par.path_JK_cov_gt
        if path_JK_cov_gt[-1] is not '/': path_JK_cov_gt+='/'
        if not os.path.exists(path_JK_cov_gt):
            os.makedirs(path_JK_cov_gt)

        # setup output path for gamma_t Jackknife covariance
        path_JK_cov_gx = self.par.path_JK_cov_gx
        if path_JK_cov_gx[-1] is not '/': path_JK_cov_gx+='/'
        if not os.path.exists(path_JK_cov_gx):
            os.makedirs(path_JK_cov_gx)

        # setup output path for gamma_t Jackknife covariance
        path_JK_cov_bf = self.par.path_JK_cov_bf
        if path_JK_cov_bf[-1] is not '/': path_JK_cov_bf+='/'
        if not os.path.exists(path_JK_cov_bf):
            os.makedirs(path_JK_cov_bf)

        # setup output path for gamma_t shot noise variance
        path_out_shot_gt = self.par.path_out_shot_gt
        if path_out_shot_gt[-1] is not '/': path_out_shot_gt+='/'
        if not os.path.exists(path_out_shot_gt):
            os.makedirs(path_out_shot_gt)

        # print feedback
        print( "Working on gamma_t calculation with bin slop=%.3f and resolution=%d:"%(self.par.bin_slop,self.par.nside) )
        print( "Running treecorr with theta=[%.1f,%.1f] over %d angular bins"%(self.par.theta_lims[0],self.par.theta_lims[1],self.par.ang_nbins) )

        # run code to get gamma_t
        for lzind in self.par.l_bins:
            # lens redshift cuts
            zl_min, zl_max = self.par.zl_bins[lzind]

            for szind in self.par.s_bins:
                # source redshift cuts
                zs_min, zs_max = self.par.zs_bins[szind]

                # give feedback on progress
                print( "  Doing: lens bin l%d [%.2f,%.2f] x source bin %d [%.2f,%.2f]"%(lzind+1,zl_min,zl_max,szind+1,zs_min,zs_max) )
                # gamma_t output directory
                gammat_out = path_out_gt+'/gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                gammax_out = path_out_gx+'/gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                extra_out = path_out_extra+'/gammat_extra_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                extra_rand_out = path_out_extra+'/gammat_rand_extra_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                randoms_gt_out = path_out_gt_rand+'/gammat_rand_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                randoms_gx_out = path_out_gx_rand+'/gammax_rand_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                boosts_out = path_out_boost+'/boost_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                cov_gammat_out = path_JK_cov_gt+'/cov_gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                cov_gammax_out = path_JK_cov_gx+'/cov_gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                cov_boosts_out = path_JK_cov_bf+'/cov_boost_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                err_gammat_out = path_JK_cov_gt+'/err_gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                shot_gammat_out = path_out_shot_gt+'/shot_noise_gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                
                # load data and setup current bin
                self.setup_run(lens_file=self.par.data_lens[lzind], lens_dir=self.par.lens_dir,
                               randoms_file=self.par.data_randoms[lzind], randoms_dir=self.par.randoms_dir, 
                               source_file=self.par.data_source[szind], 
                               lens_bin=lzind, source_bin=szind, 
                               zl_lims=[zl_min,zl_max], zs_lims=[zs_min,zs_max])

                print('Number of lenses=',len(self.ra_l))

                # random points
                if self.par.use_randoms or self.par.use_boosts:
                    print('Number of randoms=',len(self.ra_rand))
                else:
                    print('Will not use random points')
                
                # parameters to parse to treecorr
                params = [self.e1_s,self.e2_s,self.R_g,self.w_g]

                # get gamma_t for defined parameters
                (theta_res, gammat_res, gammat_total, gammat_rand, gammax_res, gammax_total, gammax_rand, 
                 cov_gammat, shot_noise_gammat, cov_boost, cov_gammax,
                 xi_im, xi_im_rand, xi_npairs, xi_npairs_rand, xi_weight, xi_weight_rand, 
                 Rg, sum_w_l, sum_w_r, 
                 boosts) = self.ggl_setup.get_gammat_and_covariance(self.ra_l, self.dec_l, ra_rand, dec_rand, self.ra_s, self.dec_s, 
                                                                    ra_rand=self.ra_rand, dec_rand=self.dec_rand, 
                                                                    params=params, low_mem=self.par.treecorr_low_mem, weights=self.weight_lens, 
                                                                    use_randoms=self.par.use_randoms, use_boosts=self.par.use_boosts)
                
                # save covariances
                #---gamma_t
                with open(cov_gammat_out,'wb') as f:
                    for line in cov_gammat.T:
                        np.savetxt(f, [line])
                gt_err = np.sqrt(np.diag(cov_gammat))
                np.savetxt(err_gammat_out, np.c_[theta_res,gt_err], header='theta, gammat_err')
                np.savetxt(shot_gammat_out, np.c_[theta_res,shot_noise_gammat], header='theta, gammat_shot_noise')
                #---boost factors
                with open(cov_boosts_out,'wb') as f:
                    for line in cov_boost.T:
                        np.savetxt(f, [line])
                #---gamma_x
                with open(cov_gammax_out,'wb') as f:
                    for line in cov_gammax.T:
                        np.savetxt(f, [line])
                
                # save results in file
                #---gamma_t
                np.savetxt(gammat_out, np.c_[theta_res,gammat_res], header='theta, gamma_t')
                np.savetxt(randoms_gt_out, np.c_[theta_res,gammat_rand], header='theta, gamma_t')
                #---gamma_x
                np.savetxt(gammax_out, np.c_[theta_res,gammax_res], header='theta, gamma_x')
                np.savetxt(randoms_gx_out, np.c_[theta_res,gammax_rand], header='theta, gamma_x')
                #---boost factors
                np.savetxt(boosts_out, np.c_[theta_res,boosts], header='theta, boost')
                #---extra stuff
                np.savetxt(extra_out, np.c_[xi_im,xi_npairs,xi_weight,Rg*np.ones(len(theta_res))], 
                                            header='xi_im, xi_npair, xi_weight, Rg')
                np.savetxt(extra_rand_out, np.c_[xi_im_rand,xi_npairs_rand,xi_weight_rand,
                                                 Rg*np.ones(len(theta_res)),sum_w_l*np.ones(len(theta_res)),sum_w_r*np.ones(len(theta_res))], 
                                                 header='xi_im, xi_npair, xi_weight, Rg, sum_w_l, sum_w_r,')

                # save results with RP subtraction and/or boost factors applied
                #---gamma_t
                if path_out_gt[-1]=='/':
                    path_out_gt_final = path_out_gt[:-1]
                else:
                    path_out_gt_final = path_out_gt
                if self.par.use_boosts:
                    path_out_gt_final += '_bf'
                if self.par.use_randoms:
                    path_out_gt_final += '_RP'
                path_out_gt_final += '/'
                gammat_total_out = path_out_gt_final+'gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                if not os.path.exists(path_out_gt_final):
                    os.makedirs(path_out_gt_final)
                np.savetxt(gammat_total_out, np.c_[theta_res,gammat_total], header='theta, gamma_t')
                #---gamma_x
                if path_out_gx[-1]=='/':
                    path_out_gx_final = path_out_gx[:-1]
                else:
                    path_out_gx_final = path_out_gx
                if self.par.use_randoms:
                    path_out_gx_final += '_RP'
                path_out_gx_final += '/'
                gammax_total_out = path_out_gx_final+'gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                if not os.path.exists(path_out_gx_final):
                    os.makedirs(path_out_gx_final)
                np.savetxt(gammax_total_out, np.c_[theta_res,gammax_total], header='theta, gamma_x')

                # give feedback on progress
                print( "  Results saved in: %s"%gammat_out )
                print( "--Done\n" )

                # clear up memory
                del self.ra_l, self.dec_l, self.ra_s, self.dec_s
                del self.e1_s,self.e2_s,self.R_g,self.w_g
                del ra_rand, dec_rand
                del self.weight_lens

        print( "Done calculating gamma_t \n" )
        return
    
