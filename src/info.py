import numpy as np
import os
from colormaps import plasma,viridis

'''
Information about paths and parameters used in pipeline. 
Define running on data (mode = data) or running or sims,
for now implemented on mice simulations (mode = mice).
Config and corresponding paths will be automatically
defined depending on the mode used. 
'''

filename_mastercat = '/global/cscratch1/sd/troxel/cats_des_y3/Y3_mastercat_v1_6_8_18_subsampled.h5'
print filename_mastercat
mode = 'data'

config_data = {
    'njk': 100,
    'bslop': 0.1,
    'nthbins': 20,
    'thlims': np.array([2.5,250.]),
    'mastercat_v': filename_mastercat[53:-3], 
    'redmagic_v': 'chi5_bright',
    'zslim_v': 'y1',
    'zs_v': 'bpz',
    'zllim_v': 'y1'
    }

config_mice = {
    'njk': 300,
    'bslop': 0.1,
    'nthbins': 20,
    'thlims': np.array([2.5,250.]),
    'version': '2',
    'redmagic_v': 'y1',
    'zslim_v': 'y1',
    'zs_v': 'bpz',
    'zllim_v': 'y1'
    }


if mode == 'data':
    config = config_data
if mode == 'data_y1sources':
    config = config_data
if mode == 'mice':
    config = config_mice

print config

paths = {}
paths['y1'] = '/Volumes/Data/y1_shear_tests/cats/jk/test_mcal_bpzmof_unblind/'
paths['runs'] =  '../runs/' 
paths['redmagic'] = '../lens_cats/redmagic/'
paths['redmagic', 'chi5_bright'] = '../lens_cats/redmagic/chi5_bright/'
paths['redmagic', 'highdens'] = '../lens_cats/redmagic/highdens/'
paths['redmagic', 'highdens_original'] = '../lens_cats/redmagic/highdens_original/'
paths['redmagic', 'y1'] = '../lens_cats/redmagic/y1/'
paths['plots'] = '../plots/'
paths['y3'] = '../../ggl_results/'
paths['y3_exp'] = '../../ggl_data/'
paths['lens'] = paths['redmagic', '%s'%config_data['redmagic_v']] + 'lens.fits'
paths['randoms'] = paths['redmagic', '%s'%config_data['redmagic_v']] + 'random.fits'
paths['mice'] = '/global/project/projectdirs/des/y3-bias/mice2/' 
paths['lens_mice'] = paths['mice'] + 'lens.fits'
paths['randoms_mice'] = paths['mice'] + 'random.fits'

paths['config_data'] = os.path.join('mastercat_%s'%config_data['mastercat_v'], 'zslim_%s'%config_data['zslim_v'], 'zs_%s'%config_data['zs_v'],
                        'redmagic_%s'%config_data['redmagic_v'], 'zllim_%s'%config_data['zllim_v'], 'njk_%d'%config_data['njk'],
                        'thbin_%0.1f_%d_%d'%(config_data['thlims'][0], config_data['thlims'][1], config_data['nthbins']),
                        'bslop_%0.1g'%config_data['bslop']) 

paths['config_mice'] = os.path.join('mice', 'v_%s'%config_mice['version'], 'zslim_%s'%config_mice['zslim_v'], 'zs_%s'%config_mice['zs_v'],
                        'redmagic_%s'%config_mice['redmagic_v'], 'zllim_%s'%config_mice['zllim_v'], 'njk_%d'%config_mice['njk'],
                        'thbin_%0.1f_%d_%d'%(config_mice['thlims'][0], config_mice['thlims'][1], config_mice['nthbins']),
                        'bslop_%0.1g'%config_mice['bslop']) 
if mode != 'data_y1sources': 
    paths['runs_config'] = os.path.join(paths['runs'], paths['config_%s'%mode]) + '/'
    paths['plots_config'] = os.path.join(paths['plots'], paths['config_%s'%mode]) + '/'
if mode == 'data_y1sources': 
    paths['runs_config'] = os.path.join(paths['runs'], paths['config_data']) + '/'
    paths['plots_config'] = os.path.join(paths['plots'], paths['config_data']) + '/'

zbins = {}
zbins['lbins'] = ['l1', 'l2', 'l3', 'l4', 'l5']
zbins['sbins'] = ['s1', 's2', 's3', 's4']
zbins['lsbins'] = [l + '_' + s for l in zbins['lbins'] for s in zbins['sbins']]
zbins['l1'] = [0.15, 0.3]
zbins['l2'] = [0.3, 0.45] 
zbins['l3'] = [0.45, 0.6] 
zbins['l4'] = [0.6, 0.75] 
zbins['l5'] = [0.75, 0.9] 
zbins['lims'] = [zbins['l1'][0], zbins['l2'][0], zbins['l3'][0], zbins['l4'][0], zbins['l5'][0], zbins['l5'][1]]
zbins['s1'] = [0.20, 0.43]
zbins['s2'] = [0.43, 0.63] 
zbins['s3'] = [0.63, 0.90] 
zbins['s4'] = [0.90, 1.30] 


plotting = {}
if mode == 'data':
    plotting['catname'] = r'Metacalibration'
if mode == 'mice':
    plotting['catname'] = r'\textsc{MICE}'

plotting['cmap'] = viridis
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

