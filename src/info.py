import numpy as np
import os
from colormaps import plasma

'''
Information about paths and parameters used in pipeline. 
'''

paths = {}

paths['runs'] =  '../runs/' 
paths['redmagic'] = '../lens_cats/redmagic/'
paths['redmagic', 'chi5_bright'] = '../lens_cats/redmagic/chi5_bright'
paths['plots'] = '../plots/'

config = {
    'njk': 100,
    'bslop': 0.1,
    'nthbins': 20,
    'thlims': np.array([2.5,250.]),
    'mastercat_v': 'v1_4_6_18', 
    'redmagic_v': 'chi5_bright',
    'zslim_v': 'y1',
    'zs_v': 'bpz',
    'zllim_v': 'y1'
    }


paths['lens'] = paths['redmagic', '%s'%config['redmagic_v']] + 'lens.fits'
paths['randoms'] = paths['redmagic', '%s'%config['redmagic_v']] + 'random.fits'

paths['runs_config'] = os.path.join(paths['runs'], 'mastercat_%s'%config['mastercat_v'], 'zslim_%s'%config['zslim_v'], 'zs_%s'%config['zs_v'],
                        'redmagic_%s'%config['redmagic_v'], 'zllim_%s'%config['zllim_v'], 'njk_%d'%config['njk'],
                        'thbin_%0.1f_%d_%d'%(config['thlims'][0], config['thlims'][1], config['nthbins']),
                        'bslop_%0.1g'%config['bslop']) + '/'

paths['plots_config'] = os.path.join(paths['plots'], 'mastercat_%s'%config['mastercat_v'], 'zslim_%s'%config['zslim_v'], 'zs_%s'%config['zs_v'],
                        'redmagic_%s'%config['redmagic_v'], 'zllim_%s'%config['zllim_v'], 'njk_%d'%config['njk'],
                        'thbin_%0.1f_%d_%d'%(config['thlims'][0], config['thlims'][1], config['nthbins']),
                        'bslop_%0.1g'%config['bslop']) + '/'

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
plotting['metacal'] = r'\textsc{Metacalibration}'
plotting['im3shape'] = r'\textsc{im3shape}'
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

