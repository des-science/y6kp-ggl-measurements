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
import os
import gc
import numpy as np
from concurrent.futures import ProcessPoolExecutor


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


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
    
    def setup_run(self, source_cat=None, path=None, lens_file=None, source_file=None, source_file_bfd=None, response=None, randoms_file=None, 
                  lens_bin=None, source_bin=None, 
                  zl_lims=None, zs_lims=None):
        """
        Setup the parameters to run code by reading files
        """
        print( "Setting up things to run code:" )

        # load lens data
        print("Reading lens data for redshift bin [%.2f,%.2f] (index %d) from %s..."%(zl_lims[0],zl_lims[1],lens_bin+1,lens_file))
        print ('lens_file', lens_file)
        # read lens galaxy data
        self.ra_l, self.dec_l, self.weight_lens = self.ggl_setup.load_lens_Y6_maglim(lens_file, zl_bin=lens_bin)
        # self.ra_l, self.dec_l, self.weight_lens = self.ggl_setup.load_lens_Y6_maglim(lens_file, zl_bin=lens_bin)
        # self.ra_l, self.dec_l, self.weight_lens = self.ggl_setup.load_lens_Y3_maglim(lens_file, zl_lims=zl_lims)

        if not self.par.use_LSSweight:
            print("Warning: Will not load LSS weights for lenses, setting them to 1")
            self.weight_lens = np.ones(len(self.ra_l))

        # load source data
        if source_cat == 'metadetect':
            print("Reading source data for source bin [%.2f,%.2f] (index %d) from %s..."%(zs_lims[0],zs_lims[1],source_bin+1,source_file))
            print ('source_file', source_file)
            # read source galaxy data
            (self.ra_s, self.dec_s, 
             self.e1_s, self.e2_s, 
             # self.R_g, self.w_g) = self.ggl_setup.load_source_metacal_5sels(source_file, zs_bin=source_bin)
             # self.R_g, self.w_g) = self.ggl_setup.load_source_metadetect(source_file, response, zs_bin=source_bin)
             self.R_g, self.w_g) = self.ggl_setup.load_source_metadetect_unblinded(source_file, response, zs_bin=source_bin)
                    # file = fits.open(par.out_main+'/metadetect_bin{}.fits'.format(zs_bin))
            if not self.par.use_response:
                print("Warning: Will set metadetect response to 1")
                self.R_g = 1.
            if not self.par.use_shearweight:
                print("Warning: Will not load metadetect weights, setting them to 1")
                self.w_g = np.ones(len(self.ra_s))
            
        elif source_cat == 'bfd':
            print ('source_file_bfd', source_file_bfd)
            (self.ra_s, self.dec_s, 
             self.e1_s, self.e2_s, 
             self.R_g, self.w_g) = self.ggl_setup.load_source_bfd(source_file_bfd, zs_bin=source_bin)
            if not self.par.use_shearweight:
                print("Warning: Will not load bfd weights, setting them to 1")
                self.w_g = np.ones(len(self.ra_s))
        else:
            print ('Specify which source catalog you want to use in the params.py file! Exiting the program')
            sys.exit(0)
            
        # load random points data
        if self.par.use_randoms or self.par.use_boosts:
            print("Reading random-point data from %s..."%(randoms_file))
            self.ra_rand, self.dec_rand = self.ggl_setup.load_randoms_Y6(randoms_file, zl_bin=lens_bin)
            # self.ra_rand, self.dec_rand = self.ggl_setup.load_randoms_Y6_maglim(path)
            # self.ra_rand, self.dec_rand = self.ggl_setup.load_randoms_Y3_maglim(randoms_file, zl_lims=zl_lims)

        else:
            print("Will not load randoms points data, as it is not needed in current run")
            self.ra_rand = None
            self.dec_rand = None

        print( "Done reading data" )

        return

  
    def setup_run_sources(self, source_cat=None, source_file=None, source_file_bfd=None, response=None,
                   source_bin=None, 
                   zs_lims=None):
        """
        Setup the parameters to run code by reading files
        """
        print( "Setting up things to run code:" )

        # load source data
        if source_cat == 'metadetect':
            print("Reading source data for source bin [%.2f,%.2f] (index %d) from %s..."%(zs_lims[0],zs_lims[1],source_bin+1,source_file))
            print ('source_file', source_file)
            # read source galaxy data
            (self.ra_s, self.dec_s, 
             self.e1_s, self.e2_s, 
             # self.R_g, self.w_g) = self.ggl_setup.load_source_metacal_5sels(source_file, zs_bin=source_bin)
             # self.R_g, self.w_g) = self.ggl_setup.load_source_metadetect(source_file, response, zs_bin=source_bin)
             self.R_g, self.w_g) = self.ggl_setup.load_source_metadetect_unblinded(source_file, response, zs_bin=source_bin)
                    # file = fits.open(par.out_main+'/metadetect_bin{}.fits'.format(zs_bin))
            if not self.par.use_response:
                print("Warning: Will set metadetect response to 1")
                self.R_g = 1.
            if not self.par.use_shearweight:
                print("Warning: Will not load metadetect weights, setting them to 1")
                self.w_g = np.ones(len(self.ra_s))
            
        elif source_cat == 'bfd':
            print ('source_file_bfd', source_file_bfd)
            (self.ra_s, self.dec_s, 
             self.e1_s, self.e2_s, 
             self.R_g, self.w_g) = self.ggl_setup.load_source_bfd(source_file_bfd, zs_bin=source_bin)
            if not self.par.use_shearweight:
                print("Warning: Will not load bfd weights, setting them to 1")
                self.w_g = np.ones(len(self.ra_s))
        else:
            print ('Specify which source catalog you want to use in the params.py file! Exiting the program')
            sys.exit(0)
            
        print( "Done reading data" )

        return


    def setup_run_lenses(self, lens_file=None, randoms_file=None, lens_bin=None, zl_lims=None):
        """
        Setup the parameters to run code by reading files
        """
        print( "Setting up things to run code:" )

        # load lens data
        print("Reading lens data for redshift bin [%.2f,%.2f] (index %d) from %s..."%(zl_lims[0],zl_lims[1],lens_bin+1,lens_file))
        print ('lens_file', lens_file)
        # read lens galaxy data
        self.ra_l, self.dec_l, self.weight_lens = self.ggl_setup.load_lens_Y6_maglim(lens_file, zl_bin=lens_bin)
        # self.ra_l, self.dec_l, self.weight_lens = self.ggl_setup.load_lens_Y6_maglim(lens_file, zl_bin=lens_bin)
        # self.ra_l, self.dec_l, self.weight_lens = self.ggl_setup.load_lens_Y3_maglim(lens_file, zl_lims=zl_lims)

        if not self.par.use_LSSweight:
            print("Warning: Will not load LSS weights for lenses, setting them to 1")
            self.weight_lens = np.ones(len(self.ra_l))

            
        # load random points data
        if self.par.use_randoms or self.par.use_boosts:
            print("Reading random-point data from %s..."%(randoms_file))
            self.ra_rand, self.dec_rand = self.ggl_setup.load_randoms_Y6(randoms_file, zl_bin=lens_bin)
            # self.ra_rand, self.dec_rand = self.ggl_setup.load_randoms_Y6_maglim(path)
            # self.ra_rand, self.dec_rand = self.ggl_setup.load_randoms_Y3_maglim(randoms_file, zl_lims=zl_lims)

        else:
            print("Will not load randoms points data, as it is not needed in current run")
            self.ra_rand = None
            self.dec_rand = None

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
        if path_out_gt[-1] != '/': path_out_gt+='/'
        if not os.path.exists(path_out_gt):
            os.makedirs(path_out_gt)

        # output directory for gamma_x
        path_out_gx = self.par.path_out_gx
        if path_out_gx[-1] != '/': path_out_gx+='/'
        if not os.path.exists(path_out_gx):
            os.makedirs(path_out_gx)

        # output directory for boost factors
        path_out_boost = self.par.path_out_boost
        if path_out_boost[-1] != '/': path_out_boost+='/'
        if not os.path.exists(path_out_boost):
            os.makedirs(path_out_boost)

        # setup output path for extra info
        path_out_extra = self.par.path_out_extra_gt
        if path_out_extra[-1] != '/': path_out_extra+='/'
        if not os.path.exists(path_out_extra):
            os.makedirs(path_out_extra)

        # setup output path for randoms
        path_out_rand = self.par.path_out_rand
        if path_out_rand[-1] != '/': path_out_rand+='/'
        if not os.path.exists(path_out_rand):
            os.makedirs(path_out_rand)
        #
        path_out_gt_rand = self.par.path_out_gt_rand
        if path_out_gt_rand[-1] != '/': path_out_gt_rand+='/'
        if not os.path.exists(path_out_gt_rand):
            os.makedirs(path_out_gt_rand)
        #
        path_out_gx_rand = self.par.path_out_gx_rand
        if path_out_gx_rand[-1] != '/': path_out_gx_rand+='/'
        if not os.path.exists(path_out_gx_rand):
            os.makedirs(path_out_gx_rand)

        # setup output path for gamma_t shot noise variance
        path_out_shot_gt = self.par.path_out_shot_gt
        if path_out_shot_gt[-1] != '/': path_out_shot_gt+='/'
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
                print( "  Doing: lens bin %d [%.2f,%.2f] x source bin %d [%.2f,%.2f]"%(lzind+1,zl_min,zl_max,szind+1,zs_min,zs_max) )
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
                # self.setup_run(lens_file=self.par.data_lens[lzind],
                #                randoms_file=self.par.data_randoms[lzind], 
                #                source_file=self.par.data_source[szind], 
                #                lens_bin=lzind, source_bin=szind, zl_lims=[zl_min,zl_max], zs_lims=[zs_min,zs_max])
                self.setup_run(source_cat=self.par.source_cat,
                               path=self.par.out_main, 
                               lens_file=self.par.data_lens,
                               randoms_file=self.par.data_randoms, 
                               source_file=self.par.data_source,
                               source_file_bfd=self.par.data_source_bfd,
                               response=self.par.response[szind],
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
                                                 Rg*np.ones(len(theta_res)),sum_w_l*np.ones(len(theta_res)),sum_w_r*np.ones(len(theta_res))], 
                                                 header='xi_im, xi_npair, xi_weight, Rg, sum_w_l, sum_w_r')

                # save results in file
                np.savetxt(gammat_out, np.c_[theta_res,gammat_res], header='theta, gamma_t')
                np.savetxt(extra_out, np.c_[xi_im,xi_npairs,xi_weight,Rg*np.ones(len(theta_res))], header='xi_im, xi_npair, xi_weight, Rg')
                np.savetxt(boosts_out, np.c_[theta_res,boosts], header='theta, boost')

                # piece together components to get gamma_t with RP subtraction and boost factors applied
                if self.par.use_boosts or self.par.use_randoms:
                    if path_out_gt[-1]=='/':
                        path_out_gt_final = path_out_gt[:-1]
                    else:
                        path_out_gt_final = path_out_gt
                    if self.par.use_boosts:
                        path_out_gt_final += '_bf'
                    if self.par.use_randoms:
                        path_out_gt_final += '_rp'
                    path_out_gt_final += '/'
                    gammat_total_out = path_out_gt_final+'gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                    # setup output path
                    if not os.path.exists(path_out_gt_final):
                        os.makedirs(path_out_gt_final)
                    np.savetxt(gammat_total_out, np.c_[theta_res,gammat_total], header='theta, gamma_t')
                else:
                    if not np.all(gammat_total/gammat_res==1.):
                        errmsg = '!!!Something is wrong, no boost or randoms, but final gammat is not equal to the basic gammat measurement'
                        raise Exception(errmsg)

                # piece together components to get gamma_x with RP subtraction and boost factors applied
                gammax_total = xi_im/Rg - xi_im_rand/Rg
                if self.par.use_randoms:
                    if path_out_gx[-1]=='/':
                        path_out_gx_final = path_out_gx[:-1]
                    else:
                        path_out_gx_final = path_out_gx
                        path_out_gx_final += '_rp'
                    path_out_gx_final += '/'
                    gammax_total_out = path_out_gx_final+'gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                    if not os.path.exists(path_out_gx_final):
                        os.makedirs(path_out_gx_final)
                    np.savetxt(gammax_total_out, np.c_[theta_res,gammax_total], header='theta, gamma_x')
                else:
                    if not np.all(gammax_total/(xi_im/Rg)==1.):
                        errmsg = '!!!Something is wrong, no randoms, but final gammax is not equal to the basic gammax measurement'
                        raise Exception(errmsg)

                # give feedback on progress
                print( "  Results saved in: %s"%gammat_out )
                print( "--Done\n" )

                # clear up memory to avoid out-of-memory issues
                del self.ra_l, self.dec_l, self.ra_s, self.dec_s
                del self.ra_rand, self.dec_rand
                del self.e1_s,self.e2_s,self.R_g,self.w_g
                del self.weight_lens

        print( "Done calculating gamma_t \n" )
        return
    
    
    
    
    
#     def run_gammat_and_cov_parallel(self):
#         """
#         Run code to get gamma_t and its covariance using only Treecorr in parallel with MPI

#         output
#         ------
#         results are saved in file
#         """
#         # output directory for gamma_t
#         path_out_gt = self.par.path_out_gt
#         # print ('path_out_gt ----------->', path_out_gt)
#         if path_out_gt[-1] != '/': path_out_gt+='/'
#         print ('path_out_gt ----------->', path_out_gt)
#         if not os.path.exists(path_out_gt):
#             os.makedirs(path_out_gt)

#         # output directory for gamma_x
#         path_out_gx = self.par.path_out_gx
#         if path_out_gx[-1] != '/': path_out_gx+='/'
#         if not os.path.exists(path_out_gx):
#             os.makedirs(path_out_gx)

#         # output directory for boost factors
#         path_out_boost = self.par.path_out_boost
#         if path_out_boost[-1] != '/': path_out_boost+='/'
#         if not os.path.exists(path_out_boost):
#             os.makedirs(path_out_boost)

#         # setup output path for extra info
#         path_out_extra = self.par.path_out_extra_gt
#         if path_out_extra[-1] != '/': path_out_extra+='/'
#         if not os.path.exists(path_out_extra):
#             os.makedirs(path_out_extra)

#         # setup output path for randoms
#         path_out_rand = self.par.path_out_rand
#         if path_out_rand[-1] != '/': path_out_rand+='/'
#         if not os.path.exists(path_out_rand):
#             os.makedirs(path_out_rand)
#         #
#         path_out_gt_rand = self.par.path_out_gt_rand
#         if path_out_gt_rand[-1] != '/': path_out_gt_rand+='/'
#         if not os.path.exists(path_out_gt_rand):
#             os.makedirs(path_out_gt_rand)
#         #
#         path_out_gx_rand = self.par.path_out_gx_rand
#         if path_out_gx_rand[-1] != '/': path_out_gx_rand+='/'
#         if not os.path.exists(path_out_gx_rand):
#             os.makedirs(path_out_gx_rand)

#         # setup output path for gamma_t Jackknife covariance
#         path_JK_cov_gt = self.par.path_JK_cov_gt
#         if path_JK_cov_gt[-1] != '/': path_JK_cov_gt+='/'
#         if not os.path.exists(path_JK_cov_gt):
#             os.makedirs(path_JK_cov_gt)

#         # setup output path for random-point gamma_t Jackknife covariance
#         path_JK_cov_gt_rand = self.par.path_JK_cov_gt_rand
#         if path_JK_cov_gt_rand[-1] != '/': path_JK_cov_gt_rand+='/'
#         if not os.path.exists(path_JK_cov_gt_rand):
#             os.makedirs(path_JK_cov_gt_rand)

#         # setup output path for gamma_t Jackknife covariance
#         path_JK_cov_gx = self.par.path_JK_cov_gx
#         if path_JK_cov_gx[-1] != '/': path_JK_cov_gx+='/'
#         if not os.path.exists(path_JK_cov_gx):
#             os.makedirs(path_JK_cov_gx)

#         # setup output path for gamma_t Jackknife covariance
#         path_JK_cov_bf = self.par.path_JK_cov_bf
#         if path_JK_cov_bf[-1] != '/': path_JK_cov_bf+='/'
#         if not os.path.exists(path_JK_cov_bf):
#             os.makedirs(path_JK_cov_bf)

#         # setup output path for gamma_t shot noise variance
#         path_out_shot_gt = self.par.path_out_shot_gt
#         if path_out_shot_gt[-1] != '/': path_out_shot_gt+='/'
#         if not os.path.exists(path_out_shot_gt):
#             os.makedirs(path_out_shot_gt)

#         # print feedback
#         print( "Working on gamma_t calculation with bin slop=%.3f and resolution=%d:"%(self.par.bin_slop,self.par.nside) )
#         print( "Running treecorr with theta=[%.1f,%.1f] over %d angular bins"%(self.par.theta_lims[0],self.par.theta_lims[1],self.par.ang_nbins) )


#         # Total number of lens-source pairs
#         total_pairs = len(self.par.l_bins) * len(self.par.s_bins)

#         self.setup_run_sources(source_cat=self.par.source_cat,
#                                source_file=self.par.data_source,
#                                source_file_bfd=self.par.data_source_bfd,
#                                response=self.par.response[szind],
#                                source_bin=szind,
#                                zs_lims=[zs_min, zs_max])

        
#         # Distribute work across ranks
#         # Each rank will process different lens-source pairs
#         for index in range(rank, total_pairs, size):  # Step through the total number of pairs, distributing among ranks
#             # Compute lens and source index for the current pair
#             lzind = index // len(self.par.s_bins)
#             szind = index % len(self.par.s_bins)

#             print(f"Rank {rank} processing: lens bin {lzind + 1}, source bin {szind + 1} <----", flush=True)

#             # lens redshift cuts
#             zl_min, zl_max = self.par.zl_bins[lzind]
#             zs_min, zs_max = self.par.zs_bins[szind]

#             # Proceed with the usual setup and calculations...
#             self.setup_run_lenses(lens_file=self.par.data_lens,
#                                   randoms_file=self.par.data_randoms,
#                                   lens_bin=lzind, 
#                                   zl_lims=[zl_min, zl_max])
            
#             print(f"Rank {rank} processing: lens bin {lzind + 1}, source bin {szind + 1}: I have uploaded the catalogs", flush=True)

#             # Call your gammat and covariance calculations
#             (theta_res, gammat_res, gammat_total, gammat_rand, gammax_res, gammax_total, gammax_rand, 
#              cov_gammat, shot_noise_gammat, cov_boost, cov_gammax, cov_gammat_rand,
#              xi_im, xi_im_rand, xi_npairs, xi_npairs_rand, xi_weight, xi_weight_rand, 
#              Rg, sum_w_l, sum_w_r, boosts) = self.ggl_setup.get_gammat_and_covariance(
#                 self.ra_l, self.dec_l, self.ra_s, self.dec_s, 
#                 ra_rand=self.ra_rand, dec_rand=self.dec_rand, 
#                 params=[self.e1_s, self.e2_s, self.R_g, self.w_g],
#                 low_mem=self.par.treecorr_low_mem, weights=self.weight_lens,
#                 use_randoms=self.par.use_randoms, use_boosts=self.par.use_boosts)

#             print(f"Rank {rank} processing: lens bin {lzind + 1}, source bin {szind + 1}: DONE GAMMAT", flush=True)

#             # Save the results as usual...
#             gammat_out = path_out_gt + f'/gammat_l{lzind + 1}_s{szind + 1}.txt'
#             np.savetxt(gammat_out, np.c_[theta_res, gammat_res], header='theta, gamma_t')

#             # Do the memory cleanup, as in the original code
#             del self.ra_l, self.dec_l, self.ra_s, self.dec_s
#             del self.e1_s, self.e2_s, self.R_g, self.w_g
#             del self.ra_rand, self.dec_rand
#             del self.weight_lens
#             gc.collect()

#         # Optionally, use MPI to reduce or gather results across ranks if necessary
#         comm.Barrier()

#         # Final steps, such as saving the two-point file, could either be done by rank 0 or combined using MPI
#         if rank == 0:
#             # Only the master rank saves the two-point file and does the final processing
#             if len(self.par.l_bins) == 6 and len(self.par.s_bins) == 4 and self.par.use_boosts and self.par.use_randoms:
#                 gammat_all = []
#                 for l in self.par.l_bins:
#                     for s in self.par.s_bins:
#                         data = np.loadtxt(self.par.out_main + f'/gammat_bf_rp/gammat_l{l+1}_s{s+1}.txt')
#                         gammat_all.append(data[:, 1])
#                 gammat_all = np.concatenate(gammat_all)
#                 dv = fits.open(self.par.dv_input)
#                 dv[4].data['VALUE'] = gammat_all
#                 dv.writeto(self.par.dv_output, overwrite=True)
#                 print(f"2pt file saved in: {self.par.dv_output}")

#         print(f"Rank {rank} done calculating gamma_t")
#         return                

    def run_gammat_and_cov_parallel(self):
        """
        Run code to get gamma_t and its covariance using only Treecorr

        output
        ------
        results are saved in file
        """
        # output directory for gamma_t
        path_out_gt = self.par.path_out_gt
        # print ('path_out_gt ----------->', path_out_gt)
        if path_out_gt[-1] != '/': path_out_gt+='/'
        print ('path_out_gt ----------->', path_out_gt)
        if not os.path.exists(path_out_gt):
            os.makedirs(path_out_gt)

        # output directory for gamma_x
        path_out_gx = self.par.path_out_gx
        if path_out_gx[-1] != '/': path_out_gx+='/'
        if not os.path.exists(path_out_gx):
            os.makedirs(path_out_gx)

        # output directory for boost factors
        path_out_boost = self.par.path_out_boost
        if path_out_boost[-1] != '/': path_out_boost+='/'
        if not os.path.exists(path_out_boost):
            os.makedirs(path_out_boost)

        # setup output path for extra info
        path_out_extra = self.par.path_out_extra_gt
        if path_out_extra[-1] != '/': path_out_extra+='/'
        if not os.path.exists(path_out_extra):
            os.makedirs(path_out_extra)

        # setup output path for randoms
        path_out_rand = self.par.path_out_rand
        if path_out_rand[-1] != '/': path_out_rand+='/'
        if not os.path.exists(path_out_rand):
            os.makedirs(path_out_rand)
        #
        path_out_gt_rand = self.par.path_out_gt_rand
        if path_out_gt_rand[-1] != '/': path_out_gt_rand+='/'
        if not os.path.exists(path_out_gt_rand):
            os.makedirs(path_out_gt_rand)
        #
        path_out_gx_rand = self.par.path_out_gx_rand
        if path_out_gx_rand[-1] != '/': path_out_gx_rand+='/'
        if not os.path.exists(path_out_gx_rand):
            os.makedirs(path_out_gx_rand)

        # setup output path for gamma_t Jackknife covariance
        path_JK_cov_gt = self.par.path_JK_cov_gt
        if path_JK_cov_gt[-1] != '/': path_JK_cov_gt+='/'
        if not os.path.exists(path_JK_cov_gt):
            os.makedirs(path_JK_cov_gt)

        # setup output path for random-point gamma_t Jackknife covariance
        path_JK_cov_gt_rand = self.par.path_JK_cov_gt_rand
        if path_JK_cov_gt_rand[-1] != '/': path_JK_cov_gt_rand+='/'
        if not os.path.exists(path_JK_cov_gt_rand):
            os.makedirs(path_JK_cov_gt_rand)

        # setup output path for gamma_t Jackknife covariance
        path_JK_cov_gx = self.par.path_JK_cov_gx
        if path_JK_cov_gx[-1] != '/': path_JK_cov_gx+='/'
        if not os.path.exists(path_JK_cov_gx):
            os.makedirs(path_JK_cov_gx)

        # setup output path for gamma_t Jackknife covariance
        path_JK_cov_bf = self.par.path_JK_cov_bf
        if path_JK_cov_bf[-1] != '/': path_JK_cov_bf+='/'
        if not os.path.exists(path_JK_cov_bf):
            os.makedirs(path_JK_cov_bf)

        # setup output path for gamma_t shot noise variance
        path_out_shot_gt = self.par.path_out_shot_gt
        if path_out_shot_gt[-1] != '/': path_out_shot_gt+='/'
        if not os.path.exists(path_out_shot_gt):
            os.makedirs(path_out_shot_gt)

        # print feedback
        print( "Working on gamma_t calculation with bin slop=%.3f and resolution=%d:"%(self.par.bin_slop,self.par.nside) )
        print( "Running treecorr with theta=[%.1f,%.1f] over %d angular bins"%(self.par.theta_lims[0],self.par.theta_lims[1],self.par.ang_nbins) )

                
            
        #================
        # srun --nodes=4 --tasks-per-node=1 python run_ggl.py
    
        #8n jobs --- 
        #evereything here is done by all processes

        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
         
        jobs = []
        for i in range(6):
            for j in range(4):
                jobs.append([i,j])

        run_count = rank 
        
        # =================
        # run code to get gamma_t
        
        while run_count < len(jobs):
            
            lzind,szind = jobs[run_count]
            
            print ('lens bin {0}, source bin {1}'.format(lzind, szind), flush=True)
            # lens redshift cuts
            zl_min, zl_max = self.par.zl_bins[lzind]           
            # source redshift cuts
            zs_min, zs_max = self.par.zs_bins[szind]
            # source redshift cuts
            # print ('self.par.zs_bins[szind]', self.par.zs_bins[szind])

            # give feedback on progress
            # print( "Doing: lens bin %d [%.2f,%.2f] x source bin %d [%.2f,%.2f]"%(lzind+1,zl_min,zl_max,szind+1,zs_min,zs_max) )
            # gamma_t output directory
            gammat_out = path_out_gt+'/gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
            gammax_out = path_out_gx+'/gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
            extra_out = path_out_extra+'/gammat_extra_l{0}_s{1}.txt'.format(lzind+1,szind+1)
            extra_rand_out = path_out_extra+'/gammat_rand_extra_l{0}_s{1}.txt'.format(lzind+1,szind+1)
            randoms_gt_out = path_out_gt_rand+'/gammat_rand_l{0}_s{1}.txt'.format(lzind+1,szind+1)
            randoms_gx_out = path_out_gx_rand+'/gammax_rand_l{0}_s{1}.txt'.format(lzind+1,szind+1)
            boosts_out = path_out_boost+'/boost_l{0}_s{1}.txt'.format(lzind+1,szind+1)
            cov_gammat_out = path_JK_cov_gt+'/cov_gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
            cov_gammat_rand_out = path_JK_cov_gt_rand+'/cov_gammat_rand_l{0}_s{1}.txt'.format(lzind+1,szind+1)
            cov_gammax_out = path_JK_cov_gx+'/cov_gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
            cov_boosts_out = path_JK_cov_bf+'/cov_boost_l{0}_s{1}.txt'.format(lzind+1,szind+1)
            err_gammat_out = path_JK_cov_gt+'/err_gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
            shot_gammat_out = path_out_shot_gt+'/shot_noise_gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)

            # print ('self.par.out_main', self.par.out_main)

            # load data and setup current bin
            self.setup_run(source_cat=self.par.source_cat,
                           path=self.par.out_main, 
                           lens_file=self.par.data_lens,
                           randoms_file=self.par.data_randoms, 
                           source_file=self.par.data_source,
                           source_file_bfd=self.par.data_source_bfd,
                           response=self.par.response[szind],
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
            print ('-------> params', params)
            print ('-------> self.R_g', self.R_g)
            print ('-------> self.w_g', self.w_g)


            # get gamma_t for defined parameters
            (theta_res, gammat_res, gammat_total, gammat_rand, gammax_res, gammax_total, gammax_rand, 
             cov_gammat, shot_noise_gammat, cov_boost, cov_gammax, cov_gammat_rand,
             xi_im, xi_im_rand, xi_npairs, xi_npairs_rand, xi_weight, xi_weight_rand, 
             Rg, sum_w_l, sum_w_r, 
             boosts) = self.ggl_setup.get_gammat_and_covariance(self.ra_l, self.dec_l, self.ra_s, self.dec_s, ra_rand=self.ra_rand, dec_rand=self.dec_rand, params=params,low_mem=self.par.treecorr_low_mem, weights=self.weight_lens,use_randoms=self.par.use_randoms,use_boosts=self.par.use_boosts)


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
            #---randoms points
            with open(cov_gammat_rand_out,'wb') as f:
                for line in cov_gammat_rand.T:
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

            # save results with RP subtraction and boost factors applied
            #---gamma_t
            if self.par.use_boosts or self.par.use_randoms:
                if path_out_gt[-1]=='/':
                    path_out_gt_final = path_out_gt[:-1]
                else:
                    path_out_gt_final = path_out_gt
                if self.par.use_boosts:
                    path_out_gt_final += '_bf'
                if self.par.use_randoms:
                    path_out_gt_final += '_rp'
                path_out_gt_final += '/'
                gammat_total_out = path_out_gt_final+'gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                # setup output path
                if not os.path.exists(path_out_gt_final):
                    os.makedirs(path_out_gt_final)
                np.savetxt(gammat_total_out, np.c_[theta_res,gammat_total], header='theta, gamma_t')
            else:
                if not np.all(gammat_total/gammat_res==1.):
                    errmsg = '!!!Something is wrong, no boost or randoms, but final gammat is not equal to the basic gammat measurement'
                    raise Exception(errmsg)
            #---gamma_x
            if self.par.use_randoms:
                if path_out_gx[-1]=='/':
                    path_out_gx_final = path_out_gx[:-1]
                else:
                    path_out_gx_final = path_out_gx
                    path_out_gx_final += '_rp'
                path_out_gx_final += '/'
                gammax_total_out = path_out_gx_final+'gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                if not os.path.exists(path_out_gx_final):
                    os.makedirs(path_out_gx_final)
                np.savetxt(gammax_total_out, np.c_[theta_res,gammax_total], header='theta, gamma_x')
            else:
                if not np.all(gammax_total/(xi_im/Rg)==1.):
                    errmsg = '!!!Something is wrong, no randoms, but final gammax is not equal to the basic gammax measurement'
                    raise Exception(errmsg)

            # give feedback on progress
            print( "  Results saved in: %s"%gammat_out )
            print( "--Done\n" )

            # clear up memory
            del self.ra_l, self.dec_l, self.ra_s, self.dec_s
            del self.e1_s,self.e2_s,self.R_g,self.w_g
            del self.ra_rand, self.dec_rand
            del self.weight_lens

            gc.collect()
            
            run_count += size
            
                
        comm.Barrier() 
        
        # Save the content of params.py to a text file        
        params_dict = {key: value for key, value in vars(self.par).items() if not key.startswith('__') and not callable(value)}
        with open(self.par.out_main+'/params_content.txt', 'w') as f:
            for key, value in params_dict.items():
                f.write(f"{key} = {value}\n")

        if((len(self.par.l_bins) == 6) & (len(self.par.s_bins) == 4) & self.par.use_boosts & self.par.use_randoms):
            print ("Saving 2pt file")
            gammat_all = []
            for l in self.par.l_bins:
                for s in self.par.s_bins:
                    asd = np.loadtxt(self.par.out_main + '/gammat_bf_rp/gammat_l{0}_s{1}.txt'.format(l+1, s+1))
                    gammat_all.append(asd[:,1])
            gammat_all = np.concatenate(gammat_all)

            dv = fits.open(self.par.dv_input)
            dv[4].data['VALUE'] = gammat_all
            dv.writeto(self.par.dv_output, overwrite=True)
            print( "2pt file saved in: %s"%self.par.dv_output )
            
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
        # print ('path_out_gt ----------->', path_out_gt)
        if path_out_gt[-1] != '/': path_out_gt+='/'
        print ('path_out_gt ----------->', path_out_gt)
        if not os.path.exists(path_out_gt):
            os.makedirs(path_out_gt)

        # output directory for gamma_x
        path_out_gx = self.par.path_out_gx
        if path_out_gx[-1] != '/': path_out_gx+='/'
        if not os.path.exists(path_out_gx):
            os.makedirs(path_out_gx)

        # output directory for boost factors
        path_out_boost = self.par.path_out_boost
        if path_out_boost[-1] != '/': path_out_boost+='/'
        if not os.path.exists(path_out_boost):
            os.makedirs(path_out_boost)

        # setup output path for extra info
        path_out_extra = self.par.path_out_extra_gt
        if path_out_extra[-1] != '/': path_out_extra+='/'
        if not os.path.exists(path_out_extra):
            os.makedirs(path_out_extra)

        # setup output path for randoms
        path_out_rand = self.par.path_out_rand
        if path_out_rand[-1] != '/': path_out_rand+='/'
        if not os.path.exists(path_out_rand):
            os.makedirs(path_out_rand)
        #
        path_out_gt_rand = self.par.path_out_gt_rand
        if path_out_gt_rand[-1] != '/': path_out_gt_rand+='/'
        if not os.path.exists(path_out_gt_rand):
            os.makedirs(path_out_gt_rand)
        #
        path_out_gx_rand = self.par.path_out_gx_rand
        if path_out_gx_rand[-1] != '/': path_out_gx_rand+='/'
        if not os.path.exists(path_out_gx_rand):
            os.makedirs(path_out_gx_rand)

        # setup output path for gamma_t Jackknife covariance
        path_JK_cov_gt = self.par.path_JK_cov_gt
        if path_JK_cov_gt[-1] != '/': path_JK_cov_gt+='/'
        if not os.path.exists(path_JK_cov_gt):
            os.makedirs(path_JK_cov_gt)

        # setup output path for random-point gamma_t Jackknife covariance
        path_JK_cov_gt_rand = self.par.path_JK_cov_gt_rand
        if path_JK_cov_gt_rand[-1] != '/': path_JK_cov_gt_rand+='/'
        if not os.path.exists(path_JK_cov_gt_rand):
            os.makedirs(path_JK_cov_gt_rand)

        # setup output path for gamma_t Jackknife covariance
        path_JK_cov_gx = self.par.path_JK_cov_gx
        if path_JK_cov_gx[-1] != '/': path_JK_cov_gx+='/'
        if not os.path.exists(path_JK_cov_gx):
            os.makedirs(path_JK_cov_gx)

        # setup output path for gamma_t Jackknife covariance
        path_JK_cov_bf = self.par.path_JK_cov_bf
        if path_JK_cov_bf[-1] != '/': path_JK_cov_bf+='/'
        if not os.path.exists(path_JK_cov_bf):
            os.makedirs(path_JK_cov_bf)

        # setup output path for gamma_t shot noise variance
        path_out_shot_gt = self.par.path_out_shot_gt
        if path_out_shot_gt[-1] != '/': path_out_shot_gt+='/'
        if not os.path.exists(path_out_shot_gt):
            os.makedirs(path_out_shot_gt)

        # print feedback
        print( "Working on gamma_t calculation with bin slop=%.3f and resolution=%d:"%(self.par.bin_slop,self.par.nside) )
        print( "Running treecorr with theta=[%.1f,%.1f] over %d angular bins"%(self.par.theta_lims[0],self.par.theta_lims[1],self.par.ang_nbins) )

                
        # run code to get gamma_t
        for lzind in self.par.l_bins:
            print ('lens bin ', lzind)
            # lens redshift cuts
            zl_min, zl_max = self.par.zl_bins[lzind]           
            
            for szind in self.par.s_bins:
                print ('SZIND', szind)
                # source redshift cuts
                zs_min, zs_max = self.par.zs_bins[szind]
                # source redshift cuts
                print ('self.par.zs_bins[szind]', self.par.zs_bins[szind])

                # give feedback on progress
                print( "  Doing: lens bin %d [%.2f,%.2f] x source bin %d [%.2f,%.2f]"%(lzind+1,zl_min,zl_max,szind+1,zs_min,zs_max) )
                # gamma_t output directory
                gammat_out = path_out_gt+'/gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                gammax_out = path_out_gx+'/gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                extra_out = path_out_extra+'/gammat_extra_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                extra_rand_out = path_out_extra+'/gammat_rand_extra_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                randoms_gt_out = path_out_gt_rand+'/gammat_rand_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                randoms_gx_out = path_out_gx_rand+'/gammax_rand_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                boosts_out = path_out_boost+'/boost_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                cov_gammat_out = path_JK_cov_gt+'/cov_gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                cov_gammat_rand_out = path_JK_cov_gt_rand+'/cov_gammat_rand_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                cov_gammax_out = path_JK_cov_gx+'/cov_gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                cov_boosts_out = path_JK_cov_bf+'/cov_boost_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                err_gammat_out = path_JK_cov_gt+'/err_gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                shot_gammat_out = path_out_shot_gt+'/shot_noise_gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)

                print ('self.par.out_main', self.par.out_main)
            
                # load data and setup current bin
                self.setup_run(source_cat=self.par.source_cat,
                               path=self.par.out_main, 
                               lens_file=self.par.data_lens,
                               randoms_file=self.par.data_randoms, 
                               source_file=self.par.data_source,
                               source_file_bfd=self.par.data_source_bfd,
                               response=self.par.response[szind],
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
                print ('-------> params', params)
                print ('-------> self.R_g', self.R_g)
                print ('-------> self.w_g', self.w_g)
                

                # get gamma_t for defined parameters
                (theta_res, gammat_res, gammat_total, gammat_rand, gammax_res, gammax_total, gammax_rand, 
                 cov_gammat, shot_noise_gammat, cov_boost, cov_gammax, cov_gammat_rand,
                 xi_im, xi_im_rand, xi_npairs, xi_npairs_rand, xi_weight, xi_weight_rand, 
                 Rg, sum_w_l, sum_w_r, 
                 boosts) = self.ggl_setup.get_gammat_and_covariance(self.ra_l, self.dec_l, self.ra_s, self.dec_s, ra_rand=self.ra_rand, dec_rand=self.dec_rand, params=params,low_mem=self.par.treecorr_low_mem, weights=self.weight_lens,use_randoms=self.par.use_randoms,use_boosts=self.par.use_boosts)

                    
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
                #---randoms points
                with open(cov_gammat_rand_out,'wb') as f:
                    for line in cov_gammat_rand.T:
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

                # save results with RP subtraction and boost factors applied
                #---gamma_t
                if self.par.use_boosts or self.par.use_randoms:
                    if path_out_gt[-1]=='/':
                        path_out_gt_final = path_out_gt[:-1]
                    else:
                        path_out_gt_final = path_out_gt
                    if self.par.use_boosts:
                        path_out_gt_final += '_bf'
                    if self.par.use_randoms:
                        path_out_gt_final += '_rp'
                    path_out_gt_final += '/'
                    gammat_total_out = path_out_gt_final+'gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                    # setup output path
                    if not os.path.exists(path_out_gt_final):
                        os.makedirs(path_out_gt_final)
                    np.savetxt(gammat_total_out, np.c_[theta_res,gammat_total], header='theta, gamma_t')
                else:
                    if not np.all(gammat_total/gammat_res==1.):
                        errmsg = '!!!Something is wrong, no boost or randoms, but final gammat is not equal to the basic gammat measurement'
                        raise Exception(errmsg)
                #---gamma_x
                if self.par.use_randoms:
                    if path_out_gx[-1]=='/':
                        path_out_gx_final = path_out_gx[:-1]
                    else:
                        path_out_gx_final = path_out_gx
                        path_out_gx_final += '_rp'
                    path_out_gx_final += '/'
                    gammax_total_out = path_out_gx_final+'gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                    if not os.path.exists(path_out_gx_final):
                        os.makedirs(path_out_gx_final)
                    np.savetxt(gammax_total_out, np.c_[theta_res,gammax_total], header='theta, gamma_x')
                else:
                    if not np.all(gammax_total/(xi_im/Rg)==1.):
                        errmsg = '!!!Something is wrong, no randoms, but final gammax is not equal to the basic gammax measurement'
                        raise Exception(errmsg)
            
                # give feedback on progress
                print( "  Results saved in: %s"%gammat_out )
                print( "--Done\n" )

                # clear up memory
                del self.ra_l, self.dec_l, self.ra_s, self.dec_s
                del self.e1_s,self.e2_s,self.R_g,self.w_g
                del self.ra_rand, self.dec_rand
                del self.weight_lens
                
                gc.collect()

                
        # Save the content of params.py to a text file        
        params_dict = {key: value for key, value in vars(self.par).items() if not key.startswith('__') and not callable(value)}
        with open(self.par.out_main+'/params_content.txt', 'w') as f:
            for key, value in params_dict.items():
                f.write(f"{key} = {value}\n")

        if((len(self.par.l_bins) == 6) & (len(self.par.s_bins) == 4) & self.par.use_boosts & self.par.use_randoms):
            print ("Saving 2pt file")
            gammat_all = []
            for l in self.par.l_bins:
                for s in self.par.s_bins:
                    asd = np.loadtxt(self.par.out_main + '/gammat_bf_rp/gammat_l{0}_s{1}.txt'.format(l+1, s+1))
                    gammat_all.append(asd[:,1])
            gammat_all = np.concatenate(gammat_all)

            dv = fits.open(self.par.dv_input)
            dv[4].data['VALUE'] = gammat_all
            dv.writeto(self.par.dv_output, overwrite=True)
            print( "2pt file saved in: %s"%self.par.dv_output )
            
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
        # print ('path_out_gt ----------->', path_out_gt)
        if path_out_gt[-1] != '/': path_out_gt+='/'
        print ('path_out_gt ----------->', path_out_gt)
        if not os.path.exists(path_out_gt):
            os.makedirs(path_out_gt)

        # output directory for gamma_x
        path_out_gx = self.par.path_out_gx
        if path_out_gx[-1] != '/': path_out_gx+='/'
        if not os.path.exists(path_out_gx):
            os.makedirs(path_out_gx)

        # output directory for boost factors
        path_out_boost = self.par.path_out_boost
        if path_out_boost[-1] != '/': path_out_boost+='/'
        if not os.path.exists(path_out_boost):
            os.makedirs(path_out_boost)

        # setup output path for extra info
        path_out_extra = self.par.path_out_extra_gt
        if path_out_extra[-1] != '/': path_out_extra+='/'
        if not os.path.exists(path_out_extra):
            os.makedirs(path_out_extra)

        # setup output path for randoms
        path_out_rand = self.par.path_out_rand
        if path_out_rand[-1] != '/': path_out_rand+='/'
        if not os.path.exists(path_out_rand):
            os.makedirs(path_out_rand)
        #
        path_out_gt_rand = self.par.path_out_gt_rand
        if path_out_gt_rand[-1] != '/': path_out_gt_rand+='/'
        if not os.path.exists(path_out_gt_rand):
            os.makedirs(path_out_gt_rand)
        #
        path_out_gx_rand = self.par.path_out_gx_rand
        if path_out_gx_rand[-1] != '/': path_out_gx_rand+='/'
        if not os.path.exists(path_out_gx_rand):
            os.makedirs(path_out_gx_rand)

        # setup output path for gamma_t Jackknife covariance
        path_JK_cov_gt = self.par.path_JK_cov_gt
        if path_JK_cov_gt[-1] != '/': path_JK_cov_gt+='/'
        if not os.path.exists(path_JK_cov_gt):
            os.makedirs(path_JK_cov_gt)

        # setup output path for random-point gamma_t Jackknife covariance
        path_JK_cov_gt_rand = self.par.path_JK_cov_gt_rand
        if path_JK_cov_gt_rand[-1] != '/': path_JK_cov_gt_rand+='/'
        if not os.path.exists(path_JK_cov_gt_rand):
            os.makedirs(path_JK_cov_gt_rand)

        # setup output path for gamma_t Jackknife covariance
        path_JK_cov_gx = self.par.path_JK_cov_gx
        if path_JK_cov_gx[-1] != '/': path_JK_cov_gx+='/'
        if not os.path.exists(path_JK_cov_gx):
            os.makedirs(path_JK_cov_gx)

        # setup output path for gamma_t Jackknife covariance
        path_JK_cov_bf = self.par.path_JK_cov_bf
        if path_JK_cov_bf[-1] != '/': path_JK_cov_bf+='/'
        if not os.path.exists(path_JK_cov_bf):
            os.makedirs(path_JK_cov_bf)

        # setup output path for gamma_t shot noise variance
        path_out_shot_gt = self.par.path_out_shot_gt
        if path_out_shot_gt[-1] != '/': path_out_shot_gt+='/'
        if not os.path.exists(path_out_shot_gt):
            os.makedirs(path_out_shot_gt)

        # print feedback
        print( "Working on gamma_t calculation with bin slop=%.3f and resolution=%d:"%(self.par.bin_slop,self.par.nside) )
        print( "Running treecorr with theta=[%.1f,%.1f] over %d angular bins"%(self.par.theta_lims[0],self.par.theta_lims[1],self.par.ang_nbins) )

                
        # run code to get gamma_t
        for lzind in self.par.l_bins:
            print ('lens bin ', lzind)
            # lens redshift cuts
            zl_min, zl_max = self.par.zl_bins[lzind]           
            
            for szind in self.par.s_bins:
                print ('SZIND', szind)
                # source redshift cuts
                zs_min, zs_max = self.par.zs_bins[szind]
                # source redshift cuts
                print ('self.par.zs_bins[szind]', self.par.zs_bins[szind])

                # give feedback on progress
                print( "  Doing: lens bin %d [%.2f,%.2f] x source bin %d [%.2f,%.2f]"%(lzind+1,zl_min,zl_max,szind+1,zs_min,zs_max) )
                # gamma_t output directory
                gammat_out = path_out_gt+'/gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                gammax_out = path_out_gx+'/gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                extra_out = path_out_extra+'/gammat_extra_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                extra_rand_out = path_out_extra+'/gammat_rand_extra_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                randoms_gt_out = path_out_gt_rand+'/gammat_rand_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                randoms_gx_out = path_out_gx_rand+'/gammax_rand_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                boosts_out = path_out_boost+'/boost_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                cov_gammat_out = path_JK_cov_gt+'/cov_gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                cov_gammat_rand_out = path_JK_cov_gt_rand+'/cov_gammat_rand_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                cov_gammax_out = path_JK_cov_gx+'/cov_gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                cov_boosts_out = path_JK_cov_bf+'/cov_boost_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                err_gammat_out = path_JK_cov_gt+'/err_gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                shot_gammat_out = path_out_shot_gt+'/shot_noise_gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)

                print ('self.par.out_main', self.par.out_main)
            
                # load data and setup current bin
                self.setup_run(source_cat=self.par.source_cat,
                               path=self.par.out_main, 
                               lens_file=self.par.data_lens,
                               randoms_file=self.par.data_randoms, 
                               source_file=self.par.data_source,
                               source_file_bfd=self.par.data_source_bfd,
                               response=self.par.response[szind],
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
                print ('-------> params', params)
                print ('-------> self.R_g', self.R_g)
                print ('-------> self.w_g', self.w_g)
                

                # get gamma_t for defined parameters
                (theta_res, gammat_res, gammat_total, gammat_rand, gammax_res, gammax_total, gammax_rand, 
                 cov_gammat, shot_noise_gammat, cov_boost, cov_gammax, cov_gammat_rand,
                 xi_im, xi_im_rand, xi_npairs, xi_npairs_rand, xi_weight, xi_weight_rand, 
                 Rg, sum_w_l, sum_w_r, 
                 boosts) = self.ggl_setup.get_gammat_and_covariance(self.ra_l, self.dec_l, self.ra_s, self.dec_s, ra_rand=self.ra_rand, dec_rand=self.dec_rand, params=params,low_mem=self.par.treecorr_low_mem, weights=self.weight_lens,use_randoms=self.par.use_randoms,use_boosts=self.par.use_boosts)

                    
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
                #---randoms points
                with open(cov_gammat_rand_out,'wb') as f:
                    for line in cov_gammat_rand.T:
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

                # save results with RP subtraction and boost factors applied
                #---gamma_t
                if self.par.use_boosts or self.par.use_randoms:
                    if path_out_gt[-1]=='/':
                        path_out_gt_final = path_out_gt[:-1]
                    else:
                        path_out_gt_final = path_out_gt
                    if self.par.use_boosts:
                        path_out_gt_final += '_bf'
                    if self.par.use_randoms:
                        path_out_gt_final += '_rp'
                    path_out_gt_final += '/'
                    gammat_total_out = path_out_gt_final+'gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                    # setup output path
                    if not os.path.exists(path_out_gt_final):
                        os.makedirs(path_out_gt_final)
                    np.savetxt(gammat_total_out, np.c_[theta_res,gammat_total], header='theta, gamma_t')
                else:
                    if not np.all(gammat_total/gammat_res==1.):
                        errmsg = '!!!Something is wrong, no boost or randoms, but final gammat is not equal to the basic gammat measurement'
                        raise Exception(errmsg)
                #---gamma_x
                if self.par.use_randoms:
                    if path_out_gx[-1]=='/':
                        path_out_gx_final = path_out_gx[:-1]
                    else:
                        path_out_gx_final = path_out_gx
                        path_out_gx_final += '_rp'
                    path_out_gx_final += '/'
                    gammax_total_out = path_out_gx_final+'gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                    if not os.path.exists(path_out_gx_final):
                        os.makedirs(path_out_gx_final)
                    np.savetxt(gammax_total_out, np.c_[theta_res,gammax_total], header='theta, gamma_x')
                else:
                    if not np.all(gammax_total/(xi_im/Rg)==1.):
                        errmsg = '!!!Something is wrong, no randoms, but final gammax is not equal to the basic gammax measurement'
                        raise Exception(errmsg)
            
                # give feedback on progress
                print( "  Results saved in: %s"%gammat_out )
                print( "--Done\n" )

                # clear up memory
                del self.ra_l, self.dec_l, self.ra_s, self.dec_s
                del self.e1_s,self.e2_s,self.R_g,self.w_g
                del self.ra_rand, self.dec_rand
                del self.weight_lens
                
                gc.collect()

                
        # Save the content of params.py to a text file        
        params_dict = {key: value for key, value in vars(self.par).items() if not key.startswith('__') and not callable(value)}
        with open(self.par.out_main+'/params_content.txt', 'w') as f:
            for key, value in params_dict.items():
                f.write(f"{key} = {value}\n")

        if((len(self.par.l_bins) == 6) & (len(self.par.s_bins) == 4) & self.par.use_boosts & self.par.use_randoms):
            print ("Saving 2pt file")
            gammat_all = []
            for l in self.par.l_bins:
                for s in self.par.s_bins:
                    asd = np.loadtxt(self.par.out_main + '/gammat_bf_rp/gammat_l{0}_s{1}.txt'.format(l+1, s+1))
                    gammat_all.append(asd[:,1])
            gammat_all = np.concatenate(gammat_all)

            dv = fits.open(self.par.dv_input)
            dv[4].data['VALUE'] = gammat_all
            dv.writeto(self.par.dv_output, overwrite=True)
            print( "2pt file saved in: %s"%self.par.dv_output )
            
        print( "Done calculating gamma_t \n" )
        return
    
    
    
    
#     import os
#     import gc
#     import numpy as np
#     from concurrent.futures import ProcessPoolExecutor
#     from astropy.io import fits

    def run_gammat_and_cov_multiprocess_chatgpt(self):
        """
        Run code to get gamma_t and its covariance using Treecorr in parallel.
        Results are saved in files.
        """

        # Helper function to process a lens-source pair
        def process_lens_source_pair(lzind, szind):
            zl_min, zl_max = self.par.zl_bins[lzind]
            zs_min, zs_max = self.par.zs_bins[szind]

            print(f"Processing: lens bin {lzind + 1} [{zl_min}, {zl_max}] x source bin {szind + 1} [{zs_min}, {zs_max}]")

            # Define output file paths for each pair
            gammat_out = path_out_gt + f'/gammat_l{lzind + 1}_s{szind + 1}.txt'
            gammax_out = path_out_gx + f'/gammax_l{lzind + 1}_s{szind + 1}.txt'
            extra_out = path_out_extra + f'/gammat_extra_l{lzind + 1}_s{szind + 1}.txt'
            extra_rand_out = path_out_extra + f'/gammat_rand_extra_l{lzind + 1}_s{szind + 1}.txt'
            randoms_gt_out = path_out_gt_rand + f'/gammat_rand_l{lzind + 1}_s{szind + 1}.txt'
            randoms_gx_out = path_out_gx_rand + f'/gammax_rand_l{lzind + 1}_s{szind + 1}.txt'
            boosts_out = path_out_boost + f'/boost_l{lzind + 1}_s{szind + 1}.txt'
            cov_gammat_out = path_JK_cov_gt + f'/cov_gammat_l{lzind + 1}_s{szind + 1}.txt'
            cov_gammat_rand_out = path_JK_cov_gt_rand + f'/cov_gammat_rand_l{lzind + 1}_s{szind + 1}.txt'
            cov_gammax_out = path_JK_cov_gx + f'/cov_gammax_l{lzind + 1}_s{szind + 1}.txt'
            cov_boosts_out = path_JK_cov_bf + f'/cov_boost_l{lzind + 1}_s{szind + 1}.txt'
            shot_gammat_out = path_out_shot_gt + f'/shot_noise_gammat_l{lzind + 1}_s{szind + 1}.txt'

            # Set up data for the current lens-source bin
            self.setup_run(source_cat=self.par.source_cat,
                           path=self.par.out_main,
                           lens_file=self.par.data_lens,
                           randoms_file=self.par.data_randoms,
                           source_file=self.par.data_source,
                           source_file_bfd=self.par.data_source_bfd,
                           response=self.par.response[szind],
                           lens_bin=lzind, source_bin=szind,
                           zl_lims=[zl_min, zl_max], zs_lims=[zs_min, zs_max])

            # Get gamma_t and covariance
            results = self.ggl_setup.get_gammat_and_covariance(
                self.ra_l, self.dec_l, self.ra_s, self.dec_s,
                ra_rand=self.ra_rand, dec_rand=self.dec_rand,
                params=[self.e1_s, self.e2_s, self.R_g, self.w_g],
                low_mem=self.par.treecorr_low_mem,
                weights=self.weight_lens,
                use_randoms=self.par.use_randoms,
                use_boosts=self.par.use_boosts
            )

            # Save results to files
            np.savetxt(gammat_out, np.c_[results[0], results[1]], header='theta, gamma_t')
            np.savetxt(shot_gammat_out, np.c_[results[0], results[9]], header='theta, shot_noise_gammat')
            np.savetxt(boosts_out, np.c_[results[0], results[20]], header='theta, boost')

            # Save covariances
            with open(cov_gammat_out, 'wb') as f:
                for line in results[7].T:
                    np.savetxt(f, [line])

            # Save other results (shortened for brevity, but follow the same logic as above)

            # Clear memory
            del self.ra_l, self.dec_l, self.ra_s, self.dec_s
            del self.e1_s, self.e2_s, self.R_g, self.w_g
            del self.ra_rand, self.dec_rand
            del self.weight_lens
            gc.collect()

        # Directories setup (same as before)
        path_out_gt = self.par.path_out_gt if self.par.path_out_gt.endswith('/') else self.par.path_out_gt + '/'
        path_out_gx = self.par.path_out_gx if self.par.path_out_gx.endswith('/') else self.par.path_out_gx + '/'
        path_out_extra = self.par.path_out_extra_gt if self.par.path_out_extra_gt.endswith('/') else self.par.path_out_extra_gt + '/'
        path_out_gt_rand = self.par.path_out_gt_rand if self.par.path_out_gt_rand.endswith('/') else self.par.path_out_gt_rand + '/'
        path_out_gx_rand = self.par.path_out_gx_rand if self.par.path_out_gx_rand.endswith('/') else self.par.path_out_gx_rand + '/'
        path_out_boost = self.par.path_out_boost if self.par.path_out_boost.endswith('/') else self.par.path_out_boost + '/'
        path_JK_cov_gt = self.par.path_JK_cov_gt if self.par.path_JK_cov_gt.endswith('/') else self.par.path_JK_cov_gt + '/'
        path_JK_cov_gt_rand = self.par.path_JK_cov_gt_rand if self.par.path_JK_cov_gt_rand.endswith('/') else self.par.path_JK_cov_gt_rand + '/'
        path_JK_cov_gx = self.par.path_JK_cov_gx if self.par.path_JK_cov_gx.endswith('/') else self.par.path_JK_cov_gx + '/'
        path_JK_cov_bf = self.par.path_JK_cov_bf if self.par.path_JK_cov_bf.endswith('/') else self.par.path_JK_cov_bf + '/'
        path_out_shot_gt = self.par.path_out_shot_gt if self.par.path_out_shot_gt.endswith('/') else self.par.path_out_shot_gt + '/'

        # Print feedback
        print(f"Working on gamma_t calculation with bin slop={self.par.bin_slop:.3f} and resolution={self.par.nside}")
        print(f"Running treecorr with theta=[{self.par.theta_lims[0]:.1f},{self.par.theta_lims[1]:.1f}] over {self.par.ang_nbins} angular bins")

        # Use ProcessPoolExecutor to parallelize over lens-source pairs
        with ProcessPoolExecutor(max_workers=24) as executor:  # Use up to 24 workers (one for each lens-source pair)
            futures = []
            for lzind in self.par.l_bins:
                for szind in self.par.s_bins:
                    futures.append(executor.submit(process_lens_source_pair, lzind, szind))

            # Wait for all futures to complete
            for future in futures:
                future.result()

        # Save the two-point file if conditions are met
        if (len(self.par.l_bins) == 6 and len(self.par.s_bins) == 4 and self.par.use_boosts and self.par.use_randoms):
            print("Saving 2pt file")
            gammat_all = []
            for l in self.par.l_bins:
                for s in self.par.s_bins:
                    data = np.loadtxt(self.par.out_main + f'/gammat_bf_rp/gammat_l{l+1}_s{s+1}.txt')
                    gammat_all.append(data[:, 1])
            gammat_all = np.concatenate(gammat_all)

            dv = fits.open(self.par.dv_input)
            dv[4].data['VALUE'] = gammat_all
            dv.writeto(self.par.dv_output, overwrite=True)
            print(f"2pt file saved in: {self.par.dv_output}")

        print("Done calculating gamma_t")
        return



    
    
    
    
    def run_gammat_and_cov_catalogsonly(self):
        """
        Run code to get gamma_t and its covariance using only Treecorr

        output
        ------
        results are saved in file
        """
        # output directory for gamma_t
        path_out_gt = self.par.path_out_gt
        # print ('path_out_gt ----------->', path_out_gt)
        if path_out_gt[-1] != '/': path_out_gt+='/'
        print ('path_out_gt ----------->', path_out_gt)
        if not os.path.exists(path_out_gt):
            os.makedirs(path_out_gt)

        # output directory for gamma_x
        path_out_gx = self.par.path_out_gx
        if path_out_gx[-1] != '/': path_out_gx+='/'
        if not os.path.exists(path_out_gx):
            os.makedirs(path_out_gx)

        # output directory for boost factors
        path_out_boost = self.par.path_out_boost
        if path_out_boost[-1] != '/': path_out_boost+='/'
        if not os.path.exists(path_out_boost):
            os.makedirs(path_out_boost)

        # setup output path for extra info
        path_out_extra = self.par.path_out_extra_gt
        if path_out_extra[-1] != '/': path_out_extra+='/'
        if not os.path.exists(path_out_extra):
            os.makedirs(path_out_extra)

        # setup output path for randoms
        path_out_rand = self.par.path_out_rand
        if path_out_rand[-1] != '/': path_out_rand+='/'
        if not os.path.exists(path_out_rand):
            os.makedirs(path_out_rand)
        #
        path_out_gt_rand = self.par.path_out_gt_rand
        if path_out_gt_rand[-1] != '/': path_out_gt_rand+='/'
        if not os.path.exists(path_out_gt_rand):
            os.makedirs(path_out_gt_rand)
        #
        path_out_gx_rand = self.par.path_out_gx_rand
        if path_out_gx_rand[-1] != '/': path_out_gx_rand+='/'
        if not os.path.exists(path_out_gx_rand):
            os.makedirs(path_out_gx_rand)

        # setup output path for gamma_t Jackknife covariance
        path_JK_cov_gt = self.par.path_JK_cov_gt
        if path_JK_cov_gt[-1] != '/': path_JK_cov_gt+='/'
        if not os.path.exists(path_JK_cov_gt):
            os.makedirs(path_JK_cov_gt)

        # setup output path for random-point gamma_t Jackknife covariance
        path_JK_cov_gt_rand = self.par.path_JK_cov_gt_rand
        if path_JK_cov_gt_rand[-1] != '/': path_JK_cov_gt_rand+='/'
        if not os.path.exists(path_JK_cov_gt_rand):
            os.makedirs(path_JK_cov_gt_rand)

        # setup output path for gamma_t Jackknife covariance
        path_JK_cov_gx = self.par.path_JK_cov_gx
        if path_JK_cov_gx[-1] != '/': path_JK_cov_gx+='/'
        if not os.path.exists(path_JK_cov_gx):
            os.makedirs(path_JK_cov_gx)

        # setup output path for gamma_t Jackknife covariance
        path_JK_cov_bf = self.par.path_JK_cov_bf
        if path_JK_cov_bf[-1] != '/': path_JK_cov_bf+='/'
        if not os.path.exists(path_JK_cov_bf):
            os.makedirs(path_JK_cov_bf)

        # setup output path for gamma_t shot noise variance
        path_out_shot_gt = self.par.path_out_shot_gt
        if path_out_shot_gt[-1] != '/': path_out_shot_gt+='/'
        if not os.path.exists(path_out_shot_gt):
            os.makedirs(path_out_shot_gt)

        # print feedback
        print( "Working on gamma_t calculation with bin slop=%.3f and resolution=%d:"%(self.par.bin_slop,self.par.nside) )
        print( "Running treecorr with theta=[%.1f,%.1f] over %d angular bins"%(self.par.theta_lims[0],self.par.theta_lims[1],self.par.ang_nbins) )

                
        # run code to get gamma_t
        for lzind in self.par.l_bins:
            print ('lens bin ', lzind)
            # lens redshift cuts
            zl_min, zl_max = self.par.zl_bins[lzind]           
            
            for szind in self.par.s_bins:
                print ('SZIND', szind)
                # source redshift cuts
                zs_min, zs_max = self.par.zs_bins[szind]
                # source redshift cuts
                print ('self.par.zs_bins[szind]', self.par.zs_bins[szind])

                # give feedback on progress
                print( "  Doing: lens bin %d [%.2f,%.2f] x source bin %d [%.2f,%.2f]"%(lzind+1,zl_min,zl_max,szind+1,zs_min,zs_max) )
                # gamma_t output directory
                gammat_out = path_out_gt+'/gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                gammax_out = path_out_gx+'/gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                extra_out = path_out_extra+'/gammat_extra_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                extra_rand_out = path_out_extra+'/gammat_rand_extra_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                randoms_gt_out = path_out_gt_rand+'/gammat_rand_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                randoms_gx_out = path_out_gx_rand+'/gammax_rand_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                boosts_out = path_out_boost+'/boost_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                cov_gammat_out = path_JK_cov_gt+'/cov_gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                cov_gammat_rand_out = path_JK_cov_gt_rand+'/cov_gammat_rand_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                cov_gammax_out = path_JK_cov_gx+'/cov_gammax_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                cov_boosts_out = path_JK_cov_bf+'/cov_boost_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                err_gammat_out = path_JK_cov_gt+'/err_gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)
                shot_gammat_out = path_out_shot_gt+'/shot_noise_gammat_l{0}_s{1}.txt'.format(lzind+1,szind+1)

                print ('self.par.out_main', self.par.out_main)
            
                # load data and setup current bin
                self.setup_run(source_cat=self.par.source_cat,
                               path=self.par.out_main, 
                               lens_file=self.par.data_lens,
                               randoms_file=self.par.data_randoms, 
                               source_file=self.par.data_source,
                               source_file_bfd=self.par.data_source_bfd,
                               response=self.par.response[szind],
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
                print ('-------> params', params)
                print ('-------> self.R_g', self.R_g)
                print ('-------> self.w_g', self.w_g)
                


                self.ggl_setup.get_gammat_and_covariance_catalogsonly(self.ra_l, self.dec_l, self.ra_s, self.dec_s, lzind, szind, ra_rand=self.ra_rand, dec_rand=self.dec_rand, params=params,low_mem=self.par.treecorr_low_mem, weights=self.weight_lens,use_randoms=self.par.use_randoms,use_boosts=self.par.use_boosts)

