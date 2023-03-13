"""
----------------------------
Main input file for Y6KP GGL
----------------------------
"""
import numpy as np
import os

config_Y3 = {
             'njk': 150,
             'bslop': 0.0,
             'nthbins': 20,
             'thlims': np.array([2.5,250.]),
             'filename_mastercat': '/global/cscratch1/sd/troxel/cats_des_y3/Y3_mastercat___UNBLIND___final_v1.0_DO_NOT_USE_FOR_2PT.h5',
             #'lens_v': 'redmagic_x40randoms_year1footprint',
             #'lens_v': 'redmagic_y1',
             #'lens_v': 'redmagic_x40randoms', # fid y3 redmagic
             'lens_v': 'maglim_x40randoms', #fid y3 maglim
             'lens_w': True,  #use LSS weights for the lenses
             #'zslim_v': 'y1',
             'zslim_v': 'som',
             #'zs_v': 'bpz',
             'zs_v': 'som',
             #'zllim_v': 'y1',
             'zllim_v': 'y3',
             'source_only_close_to_lens': True,
             'nside': 4, 
            }

def path_config(config):
    '''
    Function that defines where the measurements are saved based on the parameters defined in the config dictionary.
    '''
    return os.path.join(config['mastercat_v'],
        'zslim_%s'%config['zslim_v'],
        'zs_%s'%config['zs_v'],
        config['lens_v'],
        'zllim_%s'%config['zllim_v'],
        'lens_w_%s'%config['lens_w'],
        'njk_%d'%config['njk'],
        'thbin_%0.2f_%d_%d'%(config['thlims'][0], config['thlims'][1], config['nthbins']),
        'bslop_%0.1g'%config['bslop'],
        'source_only_close_to_lens_%s_nside%d'%(config['source_only_close_to_lens'], config['nside'])
    ) 