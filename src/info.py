import numpy as np
import os
from colormaps import plasma

'''
Information about paths and parameters used in pipeline. 
'''

root =  '/Volumes/Data/ggl_tests/'
path_cats = root + 'cats/y3/'
paths = {}

paths['runs'] = root + 'runs/' 
paths['metacal'] = path_cats + '/metacal/'
paths['metacal','v1'] = paths['metacal'] + 'v1/' 
paths['metacal','v1','cats'] = [ paths['metacal','v1'] + 'y3v02-mcal-001-combined-blind-v1-%d.fits'%n for n in range(1)]
paths['redmagic'] = path_cats + '/redmagic/'
paths['redmagic','v1'] = paths['redmagic'] + 'v1/'
paths['redmagic','v1','hidens'] = paths['redmagic','v1'] + 'y3_gold_1.0_wide_cmv02a_cmcm2_run_redmapper_v6.4.20_redmagic_highdens_0.5-10.fit'
paths['redmagic','v1','hidens_randoms'] = paths['redmagic','v1'] + 'y3_gold_1.0_wide_cmv02a_cmcm2_run_redmapper_v6.4.20_redmagic_highdens_0.5-10_randoms.fit'
paths['redmagic','v1','hidens_mask'] = paths['redmagic','v1'] + 'y3_gold_1.0_wide_cmv02a_cmcm2_run_redmapper_v6.4.20_redmagic_highdens_0.5_vlim_zmask.fit'
paths['y1'] = '/Volumes/Data/y1_shear_tests/cats/jk/test_mcal_bpzmof_unblind/'
paths['y1base'] = '/Volumes/Data/y1_shear_tests/'
paths['plots'] = root + 'plots/'
paths['nz_lens'] = paths['y1base']+'cats/redmagic/'
paths['nz_source_notomo_size'] = '/Volumes/Data/y1_shear_tests/runs/test_mcal_bpzmof_unblind/nofzs/source_size_notomo'
paths['nz_source_notomo_snr'] = '/Volumes/Data/y1_shear_tests/runs/test_mcal_bpzmof_unblind/nofzs/source_snr_notomo'
paths['sims'] = '/Volumes/Data/y1_shear_tests/runs_sims/yuedong_runs/run3/notomo/measurements/mean_measurements_l1_snotomo'
paths['cov_sims'] = '/Volumes/Data/y1_shear_tests/runs_sims/yuedong_runs/run3/notomo/covs/mean_covgt_l1_snotomo'

config = {
    'njk': 100,
    'ncpu': 1,
    'bslop': 0.1,
    'nthbins': 20,
    'thlims': np.array([2.5,250.]),
    'redmagic_v': 'y1',
    'metacal_v': 'y1',
    'zbinsource_v': 'y1'
    }


paths['runs_config'] = os.path.join(paths['runs'], 'metacal_%s'%config['metacal_v'], 'zbinsource_%s'%config['zbinsource_v'],
                        'redmagic_%s'%config['redmagic_v'], 'njk_%d'%config['njk'],
                        'thbin_%0.1f_%d_%d'%(config['thlims'][0], config['thlims'][1], config['nthbins']),
                        'bslop_%0.1g'%config['bslop']) + '/'

paths['plots_config'] = (paths['plots'] + '/metacal_%s_'%config['metacal_v'] + 'zbinsource_%s_'%config['zbinsource_v'] + 
                         'redmagic_%s_'%config['redmagic_v'] + 'njk_%d_'%config['njk']+ 
                         'thbin_%0.1f_%d_%d_'%(config['thlims'][0], config['thlims'][1], config['nthbins'])+
                         'bslop_%0.1g'%config['bslop'] + '/' )


zbins = {}
zbins['lbins'] = ['l1', 'l2', 'l3', 'l4', 'l5']
zbins['sbins'] = ['s1', 's2', 's3', 's4']
zbins['lsbins'] = [l + '_' + s for l in zbins['lbins'] for s in zbins['sbins']]
zbins['l1'] = [0.15, 0.3]
zbins['l2'] = [0.3, 0.45] 
zbins['l3'] = [0.45, 0.6] 
zbins['l4'] = [0.6, 0.75] 
zbins['l5'] = [0.75, 0.9] 
zbins['lims'] = [0.15, 0.30, 0.45, 0.60, 0.75, 0.90]
zbins['s1'] = [0.20, 0.43]
zbins['s2'] = [0.43, 0.63] 
zbins['s3'] = [0.63, 0.90] 
zbins['s4'] = [0.90, 1.30] 


plotting = {}
plotting['metacal'] = r'\textsc{Metacalibration}'
plotting['im3shape'] = r'\textsc{im3shape}'
plotting['cmap'] = plasma
plotting['redshift_l'] = [r'$0.15 < z_l < 0.30 $', r'$0.30 < z_l < 0.45$', r'$0.45 < z_l < 0.60$', r'$0.60 < z_l < 0.75 $', r'$0.75 < z_l < 0.90 $']
plotting['redshift_s'] = [r'$0.20 < z_s < 0.43$', r'$0.43 < z_s < 0.63  $', r'$0.63 < z_s < 0.90  $', r'$0.90 < z_s < 1.30  $']
plotting['titles_redmagic'] = ['redMaGiC HiDens', 'redMaGiC HiDens', 'redMaGiC HiDens', 'redMaGiC HiLum', 'redMaGiC HigherLum']
plotting['th_limit'] = [64.,40.,30., 24., 21.]


source_nofz_pars = {}
source_nofz_pars['dzs','size'] = [-0.00071855,0.003097] # true - bpz for [low, high] size
source_nofz_pars['dzs','snr'] = [0.032083,-0.003159] # true - bpz for [low, high] snr
source_nofz_pars['dzs_sigma'] = 1.6 #sqrt(2) times the non-tomographic uncertainty as estimated in Hoyle et al. using COSMOS
source_nofz_pars['thetamin'] = 64.

sysmaps = {}
sysmaps['nside'] = 4096
sysmaps['nested_bool'] = False # means that it is ring
sysmaps['airmass'] = 'AIRMASS_coaddweights3_mean'
sysmaps['count'] = 'count__fracdet'
sysmaps['exptime'] = 'EXPTIME__total'
sysmaps['fwhm'] = 'FWHM_MEAN_coaddweights3_mean'
sysmaps['maglimit'] = 'maglimit3__'
sysmaps['skybrite'] = 'SKYBRITE_coaddweights3_mean'
sysmaps['skysigma'] = 'SKYSIGMA_coaddweights3_mean'

