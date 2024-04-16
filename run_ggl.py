import sys
import os
import numpy as np
import yaml
import pandas as pd
from astropy.io import fits
import gc
from mpi4py import MPI
import setup
import ggl


def printlog(msg, file_name):
    print(msg)
    with open(file_name, 'a+') as f:
        f.write(f'{msg}\n')


class run_GGL(object):


    def __init__(self, config=None):
        
        self.par = config
        return
    
    
    def create_output_dir(self):
        dir_out = self.par['dir_out']
        print(f'Saving GGL measurements in {dir_out}')
        if not os.path.exists(dir_out):
            os.makedirs(dir_out)
        if self.par['compute_covariance']:
            if not os.path.exists(f'{dir_out}/covariance/'):
                os.makedirs(f'{dir_out}/covariance/')
        return
    

    def save_results_2pt_file(self):
        
        dir_out = self.par['dir_out']
        
        theta = []
        gammat = []
        for l_zbin in self.par['l_bins']:
            for s_zbin in self.par['s_bins']:
                ggl_results = np.genfromtxt(f'{dir_out}/ggl_l{l_zbin}_s{s_zbin}.txt', names=True)
                theta.append(ggl_results['theta'])
                gammat.append(ggl_results['gammat_bf_rp'])
        theta = np.concatenate(theta)                
        gammat = np.concatenate(gammat)

        dv_out = f'{dir_out}/'+self.par['dv_output']
        with fits.open(self.par['dv_input']) as dv:
            dv[4].data['ANG'] = theta
            dv[4].data['VALUE'] = gammat
            dv.writeto(dv_out, overwrite=True)
        
        print(f'2pt file saved in {dv_out}')
        return
    
    
    def save_results_pkl_file(self):
        
        dir_out = self.par['dir_out']
        
        ggl_file0 = f'{dir_out}/ggl_l'+str(self.par['l_bins'][0])+'_s'+str(self.par['s_bins'][0])+'.txt'
        with open(ggl_file0, 'r') as f:
            columns = f.readline().split()[1:]
        ggl_data = pd.DataFrame(columns=columns)

        l_zbin_col = []
        s_zbin_col = []

        for l_zbin in self.par['l_bins']:
            for s_zbin in self.par['s_bins']:

                l_zbin_col.append([l_zbin]*self.par['ang_nbins'])
                s_zbin_col.append([s_zbin]*self.par['ang_nbins'])

                ggl_file = f'{dir_out}/ggl_l{l_zbin}_s{s_zbin}.txt'
                ggl_data_bin = pd.DataFrame(np.genfromtxt(ggl_file, names=True))

                ggl_data = pd.concat([ggl_data, ggl_data_bin])

        ggl_data['l_zbin'] = np.concatenate(l_zbin_col)
        ggl_data['s_zbin'] = np.concatenate(s_zbin_col)
        ggl_data.reset_index()
        
        cov_types = ['gammat', 'gammax']
        if self.par['use_randoms']: cov_types.append('rand')
        if self.par['use_boost']: cov_types.append('boost')
        
        for cov_type in cov_types:
            err = []
            for l_zbin in self.par['l_bins']:
                for s_zbin in self.par['s_bins']:
                    cov_file = f'{dir_out}/covariance/cov_l{l_zbin}_s{s_zbin}_{cov_type}.txt'
                    cov = np.loadtxt(cov_file)
                    err.append(np.sqrt(np.diag(cov.T)))
        
            ggl_data[f'err_{cov_type}'] = np.concatenate(err)
        
        
        pickle_out = f'{dir_out}/'+self.par['dv_output'][:-5]+'.pkl'
        ggl_data.to_pickle(pickle_out)
        
        print(f'pkl file saved in {pickle_out}')
        return
    
    
    def save_params(self):
        
        params_dict = {key: value for key, value in self.par.items()
            if not key.startswith('__') and not callable(value)}
        
        params_file_out = self.par['dir_out']+'/params_content.txt'
        with open(params_file_out, 'w') as f:
            for key, value in params_dict.items():
                f.write(f'{key} = {value}\n')
                
        print(f'Parameters saved in {params_file_out}')
        return
    
    
    def printlog_bins_info(self, l_zbin, s_zbin):
        
        bins_file = self.par['dir_out']+f'/info_l{l_zbin}_s{s_zbin}.txt'
        if os.path.exists(bins_file):
            os.remove(bins_file)

        printlog('Num lenses = {:d}'.format(len(self.ra_l)), bins_file)
        printlog('Num sources = {:d}'.format(len(self.ra_s)), bins_file)
        printlog('Num randoms = {:d}'.format(len(self.ra_r)), bins_file)
        printlog('Sum weights sources = {:f}'.format(np.sum(self.w_s)), bins_file)
        printlog('Sum weights lenses = {:f}'.format(np.sum(self.w_l)), bins_file)
        
        if (self.par['source_cat'] == 'bfd') & (self.par['ggl_bfd_approximate'] == False):
            printlog('Average Q1 = {:f}'.format(np.average(self.Q0, weights=self.w_s)), bins_file)
            printlog('Average Q2 = {:f}'.format(np.average(self.Q1, weights=self.w_s)), bins_file)
        else:
            printlog('Response = {:f}'.format(self.R), bins_file)
            printlog('Average e1 = {:f}'.format(np.average(self.e1, weights=self.w_s)), bins_file)
            printlog('Average e2 = {:f}'.format(np.average(self.e2, weights=self.w_s)), bins_file)

        return
       

    def load_catalogs(self, l_zbin, s_zbin, zl_lims, zs_lims):
        """Load lens-sources bins pair data"""

        # load lens data
        lens_cat = self.par['lens_cat']
        print('Reading {} lens data, redshift bin {:d} [{:.3f},{:.3f}] from {}'.format(
            lens_cat, l_zbin+1, zl_lims[0], zl_lims[1], self.par[f'data_lens_{lens_cat}']))
        
        if lens_cat == 'maglim':
            (self.ra_l, self.dec_l, self.w_l) = setup.load_lens_Y6_maglim(self.par['data_lens_maglim'], self.par['data_LSSweights'], l_zbin)
        elif lens_cat == 'maglim_y3':
            (self.ra_l, self.dec_l, self.w_l) = setup.load_lens_Y3_maglim(self.par['data_lens_maglim_y3'], zl_lims)
        else:
            print('Specify which lens catalog you want to use in the configuration file! Exiting the program')
            sys.exit(0)

        # load source data
        source_cat = self.par['source_cat']
        print('Reading {} source data, redshift bin {:d} [{:.3f},{:.3f}] from {}'.format(
            source_cat, s_zbin+1, zs_lims[0], zs_lims[1], self.par[f'data_source_{source_cat}']))
            
        if self.par['source_cat'] == 'metadetect':
            (self.ra_s, self.dec_s, self.e1, self.e2, self.R, self.w_s) = setup.load_source_metadetect(self.par['data_source_metadetect'], s_zbin)
        
        elif self.par['source_cat'] == 'bfd':
            (self.ra_s, self.dec_s, self.e1, self.e2, 
             self.P, self.Q0, self.Q1, 
            self.R00, self.R01, self.R11, self.w_s) = setup.load_source_bfd(self.par['data_source_bfd'],self.par['data_source_bfd_binning'], 
                                                                            self.par['data_source_bfd_mask'], s_zbin)
            self.R = 1.0
            
        elif self.par['source_cat'] == 'metacal':
            (self.ra_s, self.dec_s, self.e1, self.e2, self.R, self.w_s) = setup.load_source_metacal(self.par['data_source_metacal'], s_zbin)

        else:
            print('Specify which source catalog you want to use in the params.py file! Exiting the program')
            sys.exit(0)

        # load random points data
        if (self.par['use_randoms'] or self.par['use_boost'] or self.par['compute_covariance']):
            print('Reading {} random-point data, redshift bin {:d} from {}'.format(lens_cat, l_zbin+1, self.par[f'data_randoms_{lens_cat}']))

            if lens_cat == 'maglim':
                (self.ra_r, self.dec_r) = setup.load_randoms_Y6(self.par['data_randoms_maglim'], l_zbin)
            elif lens_cat == 'maglim_y3':
                (self.ra_r, self.dec_r) = setup.load_randoms_Y3(self.par['data_randoms_maglim_y3'], zl_lims)

        else:
            print('Will not load randoms points data, as it is not needed in current run')
            self.ra_r = None
            self.dec_r = None

        if not self.par['use_LSSweight']:
            print('WARNING: Will not use LSS weights for lenses, setting them to 1')
            self.w_l = np.ones(len(self.ra_l))

        if not self.par['use_response']:
            print('WARNING: Will set Response to 1')
            self.R = 1.0

        if not self.par['use_shearweight']:
            print('WARNING: Will not use source weights, setting them to 1')
            self.w_s = np.ones(len(self.ra_s))

        return


    def run_ggl(self, l_zbin, s_zbin):
        """Run GGL measurements for a source-lens bins pair"""

        # lens redshift cuts
        zl_min, zl_max = self.par['zl_bins'][l_zbin]
        # source redshift cuts
        zs_min, zs_max = self.par['zs_bins'][s_zbin]

        print('\nLens bin {:d} [{:.3f},{:.3f}] x source bin {:d} [{:.3f},{:.3f}]'.format(
            l_zbin+1, zl_min, zl_max, s_zbin+1, zs_min, zs_max))
        sys.stdout.flush()

        # load data for current bins pair
        self.load_catalogs(l_zbin, s_zbin, zl_lims=[zl_min, zl_max], zs_lims=[zs_min, zs_max])
        self.printlog_bins_info(l_zbin, s_zbin)

        dir_out = self.par['dir_out']
        ggl_file_out = f'{dir_out}/ggl_l{l_zbin}_s{s_zbin}.txt'
        cov_file_out = f'{dir_out}/covariance/cov_l{l_zbin}_s{s_zbin}.txt'
        
        if (self.par['source_cat'] == 'bfd') & (self.par['ggl_bfd_approximate'] == False):
            
            # run BFD ggl measurements for current bins pair
            ggl.get_ggl_bfd(self.ra_l, self.dec_l, self.w_l, self.ra_s, self.dec_s, self.w_s, 
                               self.Q0, self.Q1, self.R00, self.R01, self.R11,
                               self.ra_r, self.dec_r, units='deg', 
                               theta_lims=self.par['theta_lims'], nbins=self.par['ang_nbins'], sep_units='arcmin',
                               bin_slop=self.par['bin_slop'], low_mem=self.par['treecorr_low_mem'],
                               use_randoms=self.par['use_randoms'], use_boost=self.par['use_boost'],
                               compute_cov=self.par['compute_covariance'], npatch=self.par['n_jck'], 
                               ggl_file_out=ggl_file_out, cov_file_out=cov_file_out)
        
        else:
            # run ggl measurements for current bins pair
            ggl.get_ggl(self.ra_l, self.dec_l, self.w_l, self.ra_s, self.dec_s, self.w_s, 
                           self.e1, self.e2, self.R,
                           self.ra_r, self.dec_r, units='deg',
                           theta_lims=self.par['theta_lims'], nbins=self.par['ang_nbins'], sep_units='arcmin',
                           bin_slop=self.par['bin_slop'], low_mem=self.par['treecorr_low_mem'],
                           use_randoms=self.par['use_randoms'], use_boost=self.par['use_boost'],
                           compute_cov=self.par['compute_covariance'], npatch=self.par['n_jck'], 
                           ggl_file_out=ggl_file_out, cov_file_out=cov_file_out)
            
        # clear up memory
        del self.ra_l, self.dec_l, self.w_l, self.ra_s, self.dec_s, self.w_s, self.ra_r, self.dec_r
        if (self.par['source_cat'] == 'bfd') & (self.par['ggl_bfd_approximate'] == False):
            del self.Q0, self.Q1, self.R00, self.R01, self.R11
        else:
            del self.e1, self.e2, self.R

        gc.collect()
        return

    
    def run_ggl_measurements(self):
        """Run GGL measurement"""
        
        print('\nComputing GGL with bin slop={:.3f}'.format(self.par['bin_slop']))
        print('Running TreeCorr with theta=[{:.1f},{:.1f}] over {:d} angular bins'.format(
                self.par['theta_lims'][0], self.par['theta_lims'][1], self.par['ang_nbins']))
        
        self.create_output_dir()

        if self.par['run_parallel']:

            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()

            # Define the total number of lens and source bins
            num_lens_bins = len(self.par['l_bins'])
            num_source_bins = len(self.par['s_bins'])

            # Calculate the chunks per process
            chunks_per_process = (num_lens_bins * num_source_bins) // size

            # Calculate the start and end chunk for this process
            start_chunk = rank * chunks_per_process
            end_chunk = start_chunk + chunks_per_process

            # Iterate over the assigned chunks
            for chunk in range(start_chunk, end_chunk):
                l_zbin = chunk // num_source_bins
                s_zbin = chunk % num_source_bins
                
                self.run_ggl(l_zbin, s_zbin)
        else:
            
            for l_zbin in self.par['l_bins']:
                for s_zbin in self.par['s_bins']:
    
                    self.run_ggl(l_zbin, s_zbin)
                        
        print('\nDone running GGL measurements')
        
        # Save the content of params.py to a text file
        self.save_params()
        
        # gather all results in a single file
        self.save_results_pkl_file()
        
        if ((len(self.par['l_bins']) == 6) & (len(self.par['s_bins']) == 4) & self.par['use_boost'] & self.par['use_randoms']):
            self.save_results_2pt_file()
           
        print('\nWARNING: all the results saved are UNBLINDED, need to run blinding (see README)\n')

        return
    


if __name__ == "__main__":

    # Load configuration file with the parameters for the GGL run
    if len(sys.argv) > 1:
        config_file = sys.argv[1]
    else:
        config_file = './input/params.yaml'
        print(f'\nWARNING: No custom configuration file provided, loading parameters from {config_file}')
    config = yaml.safe_load(open(config_file, 'r'))
    
    run = run_GGL(config)
    run.run_ggl_measurements()
    
    exit()