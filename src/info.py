import numpy as np
import os
from colormaps import plasma,viridis

'''
Information about paths and parameters used in pipeline. 
'''

"""
BASIC
---------------------------------
Define the basic parameters here. Define running on data (mode = 'data') 
or running or sims, for now implemented on mice simulations (mode = 'mice').
Config and corresponding paths will be automatically defined depending on the mode used. 

BLINDING Instructions: 
If you are running measurements on data, you need to blind them. Steps:

  1) Set the following below:
     blind = 1
     run = 1
     savetwopoint = 1
     plot = 0

  2) Run the scripts normally (i.e. run call_ggl_y3kp). 
  This will save a TwoPointFile with the unblinded measurements. 

  3) Source the setup file of cosmosis (i.e. source cosmosis/config/setup-cosmosis-nersc, 
  using your own path to the file). This will enable the cosmosis environment. 

  4) Run the blind.py script (you might need to pull from 2pt_pipeline repository first, 
  and edit the path to that directory in blind.py). The blind.py script will blind the 
  measurements using cosmosis, will save a twopoint file with the blinded measurements 
  and will delete the unblinded ones. 

  5) Then, to plot the blinded measurements, set 
     blind = 1
     run =  0
     savetwopoint = 0
     plot = 1
     and run call_ggl_y3kp.py again. 
     
Then you are done. All this process is only need for the gammat measurements, not the systematics tests. 
"""

basic = {
    'mode':'data',
    'blind': 1,
    'run': 0,
    'savetwopoint': 1,
    'plot': 0,
    'pool': 1,
    'computer': 'nersc'  #can be 'nersc', 'local', 'midway', etc
}

if basic['pool']:
    basic['num_threads'] = 1
    basic['Ncores'] = 7
    
if not basic['pool']: basic['num_threads'] = 10


"""
CONFIGURATION
---------------------------------
Define the parameters of your run here. 
If you are running on data edit config_data, 
and if you are using MiCE simulations use
config_mice.
"""

config_data = {
    'njk': 150,
    'bslop': 0.1,
    'nthbins': 20,
    'thlims': np.array([2.5,250.]),
    'filename_mastercat': '/project/projectdirs/des/www/y3_cats/Y3_mastercat_12_3_19.h5',
    'lens_v': 'redmagic',
    #'lens_v': 'maglim',
    'lens_w': True,  #use LSS weights for the lenses
    'zslim_v': 'som',
    'zs_v': 'bpz',
    'zllim_v': 'y3',
    }

config_data['mastercat_v'] = config_data['filename_mastercat'][-23:-3]
if basic['computer']=='midway':
    config_data['filename_mastercat'] = '/project2/chihway/data/des_y3_catalogs/y3kp_sample/' + config_data['mastercat_v'] + '.h5'

print config_data['filename_mastercat']
config_mice = {
    'njk': 300,
    'bslop': 0.1,
    'nthbins': 20,
    'thlims': np.array([2.5,250.]),
    'version': '2',
    'lens_v': 'y1',
    'zslim_v': 'y1',
    'zs_v': 'bpz',
    'zllim_v': 'y1'
    }

if basic['mode'] == 'data':
    config = config_data
if basic['mode'] == 'mice':
    config = config_mice

print '\nChosen configuration:\n--------------------------\n', config

"""
PATHS
---------------------------------
Define the paths dictionary. This dictionary is imported
in the other scripts. Add more paths as necessary.
"""

paths = {}
paths['runs'] =  '../runs/' 
paths['plots'] = '../plots/'
paths['yaml'] = 'cats.yaml' 
paths['y3'] = '../../ggl_results/'
paths['y3_exp'] = '../../ggl_data/'
if basic['computer'] == 'nersc':
    paths['y3_sysmap'] = '/global/project/projectdirs/des/ggl/systematic_maps/'
    paths['lens_cats'] = '/global/project/projectdirs/des/ggl/lens_cats/'

if basic['computer'] == 'midway':
    paths['lens_cats'] = '../cats/'

if 'redmagic' in config['lens_v']:
    paths['lens_nz'] = '../simulated_dvs/sim_fiducial_redmagic_sompzv0.132_covDec2019.fits'
if 'maglim' in config['lens_v']:
    paths['lens_nz'] = '../simulated_dvs/sim_fiducial_maglim_sompzv0.132_covDec2019.fits'

paths['source_nz'] = paths['lens_nz']
#paths['source_nz'] = 'y1_2pt_NG_mcal_1110.fits'
print paths['lens_nz']

if basic['mode'] == 'data':
    paths['lenscats'] = paths['lens_cats'] + '%s/%s/njk_%d/'%(config['mastercat_v'], config['lens_v'],config['njk']) 
    paths['lens'] = paths['lenscats'] + 'lens.fits'
    paths['randoms'] = paths['lenscats'] + 'random.fits'
    print '--------------------------\nUsing lens file in:\n', paths['lens'] 
    print '--------------------------\nUsing randoms file in:\n', paths['randoms'] 

    paths['config_data'] = os.path.join('mastercat_%s'%config_data['mastercat_v'], 'zslim_%s'%config_data['zslim_v'], 'zs_%s'%config_data['zs_v'],
                                        config_data['lens_v'], 'zllim_%s'%config_data['zllim_v'], 'lens_w_%s'%config_data['lens_w'], 'njk_%d'%config_data['njk'],
                                        'thbin_%0.1f_%d_%d'%(config_data['thlims'][0], config_data['thlims'][1], config_data['nthbins']),
                                        'bslop_%0.1g'%config_data['bslop']) 

    paths['theory_size_all_covmat'] = paths['runs']+paths['config_data'] + '/size/theory_size_all_covmat.fits'
    paths['hist_n_of_z_lenses_witherr_size'] = paths['runs']+paths['config_data'] +'/size/hist_n_of_z_lenses_witherr_size'
    paths['hist_n_of_z_low_size'] = paths['runs']+paths['config_data']+'/size/hist_n_of_z_low_size'
    paths['hist_n_of_z_high_size'] = paths['runs']+paths['config_data']+'/size/hist_n_of_z_high_size'
    
if basic['mode'] == 'mice':
    paths['mice'] = '/global/project/projectdirs/des/y3-bias/mice2/' 
    paths['lens_mice'] = paths['mice'] + 'lens.fits'
    paths['randoms_mice'] = paths['mice'] + 'random.fits'
    paths['config_mice'] = os.path.join('mice', 'v_%s'%config_mice['version'], 'zslim_%s'%config_mice['zslim_v'], 'zs_%s'%config_mice['zs_v'],
                                        'redmagic_%s'%config_mice['lens_v'], 'zllim_%s'%config_mice['zllim_v'], 'njk_%d'%config_mice['njk'],
                                        'thbin_%0.1f_%d_%d'%(config_mice['thlims'][0], config_mice['thlims'][1], config_mice['nthbins']),
                                        'bslop_%0.1g'%config_mice['bslop']) 

# Where we save the runs and plots for one particular configuration:
paths['runs_config'] = os.path.join(paths['runs'], paths['config_%s'%basic['mode']]) + '/'
paths['plots_config'] = os.path.join(paths['plots'], paths['config_%s'%basic['mode']]) + '/'


print paths['runs_config']
"""
ZBINS
---------------------------------
Define the zbins dictionary. This dictionary is imported
in the other scripts and it defines the number of lens and source
redshift bins and their limits. 
"""

zbins = {}
if 'redmagic' in config['lens_v'] or basic['mode']=='mice':
    zbins['lbins'] = ['l1', 'l2', 'l3', 'l4', 'l5']
if 'maglim' in config['lens_v']:
    zbins['lbins'] = ['l1', 'l2', 'l3', 'l4', 'l5', 'l6']
#zbins['sbins'] = ['s2', 's3', 's4'] 
zbins['sbins'] = ['s1', 's2', 's3', 's4']
#zbins['sbins'] = ['s3']
zbins['lsbins'] = [l + '_' + s for l in zbins['lbins'] for s in zbins['sbins']]
zbins['sys'] = [0.2, 1.20]

if config['zllim_v'] == 'y1':
    zbins['l1'] = [0.15, 0.30]
    zbins['l2'] = [0.30, 0.45] 
    zbins['l3'] = [0.45, 0.60] 
    zbins['l4'] = [0.60, 0.75] 
    zbins['l5'] = [0.75, 0.90] 
    zbins['lims'] = [zbins['l1'][0], zbins['l2'][0], zbins['l3'][0], zbins['l4'][0], zbins['l5'][0], zbins['l5'][1]]

if config['zllim_v'] == 'y3':
    zbins['l1'] = [0.15, 0.35]
    zbins['l2'] = [0.35, 0.50] 
    zbins['l3'] = [0.50, 0.65] 
    zbins['l4'] = [0.65, 0.80] 
    zbins['l5'] = [0.80, 0.90] 
    zbins['lims'] = [zbins['l1'][0], zbins['l2'][0], zbins['l3'][0], zbins['l4'][0], zbins['l5'][0], zbins['l5'][1]]

if 'maglim' in config['lens_v']:
    zbins['l1'] = [0.20, 0.40]
    zbins['l2'] = [0.40, 0.55] 
    zbins['l3'] = [0.55, 0.70] 
    zbins['l4'] = [0.70, 0.85] 
    zbins['l5'] = [0.85, 0.95] 
    zbins['l6'] = [0.95, 1.05] 
    zbins['lims'] = [zbins['l1'][0], zbins['l2'][0], zbins['l3'][0], zbins['l4'][0], zbins['l5'][0], zbins['l6'][0],zbins['l6'][1]]


if config['zslim_v'] == 'som':
    # SOM bins don't have explicit limits
    zbins['s1'] = [0, 0]
    zbins['s2'] = [1, 1] 
    zbins['s3'] = [2, 2] 
    zbins['s4'] = [3, 3]
if config['zslim_v'] == 'y1':
    zbins['s1'] = [0.20, 0.43]
    zbins['s2'] = [0.43, 0.63] 
    zbins['s3'] = [0.63, 0.90] 
    zbins['s4'] = [0.90, 1.30]
    zbins['source_lims'] = [zbins['s1'][0], zbins['s2'][0], zbins['s3'][0], zbins['s4'][0], zbins['s4'][1]]

"""
PLOTTING
---------------------------------
Define the plotting dictionary. This dictionary is imported
in the other scripts and it defines useful quantites that are used
accross several plots, to ensure consistency. 
"""

plotting = {}
if basic['mode'] == 'data':
    plotting['catname'] = r'Metacalibration ' + config['mastercat_v'][0:2] + ' ' + config['lens_v']
if basic['mode'] == 'mice':
    plotting['catname'] = r'\textsc{MICE}'

plotting['latex'] = False
plotting['cmap'] = viridis
plotting['redshift_l'] = [r'$%0.2f < z_l < %0.2f $'%(zbins['lims'][i], zbins['lims'][i+1]) for i in range(len(zbins['lims'])-1)]
#plotting['th_limit'] = [64.,40.,30., 24., 21.] # 12 Mpc/h 
#plotting['th_limit'] = [42.67, 26.67 ,20., 16., 14.] # 8 Mpc/h 
plotting['th_limit'] = [21.33, 13.33 , 10., 8., 7., 6.] # 4 Mpc/h #check values for last bin
if config['zslim_v'] == 'y1':
    plotting['redshift_s'] = [r'$%0.2f < z_s < %0.2f $'%(zbins['source_lims'][i], zbins['source_lims'][i+1]) for i in range(len(zbins['source_lims'])-1)]
if config['zslim_v'] == 'som':
    plotting['redshift_s'] = ['Bin 1', 'Bin 2', 'Bin 3', 'Bin 4']

plotting['titles_redmagic'] = ['redMaGiC HiDens', 'redMaGiC HiDens', 'redMaGiC HiDens', 'redMaGiC HiLum', 'redMaGiC HiLum']


"""
SIZE AND S/N NOISE TESTS
---------------------------------
Defines a dictionary used for the size and S/N. 
Currently values from y1. Needs to be updated for Y3.
"""

source_nofz_pars = {}
source_nofz_pars['dzs','size'] = [-0.00071855,0.003097] # true - bpz for [low, high] size
source_nofz_pars['dzs','snr'] = [0.032083,-0.003159] # true - bpz for [low, high] snr
source_nofz_pars['dzs_sigma'] = 1.6 #sqrt(2) times the non-tomographic uncertainty as estimated in Hoyle et al. using COSMOS
source_nofz_pars['thetamin'] = 64.

"""
SYSTEMATICS MAPS TESTS
---------------------------------
Defines a dictionary used for the systematics maps tests.
Currently names for y1. 
"""

sysmaps = {}
sysmaps['nside'] = 4096
sysmaps['nested_bool'] = True # means that it is ring
# sysmaps['airmass'] = 'AIRMASS_coaddweights3_mean'

sysmaps['airmass'] = 'AIRMASS.WMEAN_EQU'

sysmaps['count'] = 'count__fracdet'
sysmaps['exptime'] = 'EXPTIME__total'
sysmaps['fwhm'] = 'FWHM.WMEAN_EQU'
sysmaps['maglimit'] = 'maglim_EQU'
sysmaps['skybrite'] = 'SKYBRITE.WMEAN_EQU'
sysmaps['skysigma'] = 'SKYSIGMA_coaddweights3_mean'
sysmaps['separate_jk'] = False
sysmaps['separate_jk_njk'] = 75

if sysmaps['separate_jk']:
    config['njk'] = sysmaps['separate_jk_njk']
