import os
from info import paths, config

path_2pt_pipeline = '../../2pt_pipeline/'


def blind_measurements(name):
        """
        Run blinding script from 2pt pipeline repository. 
        In order to run this, cosmosis environment needs to be sourced.
        It needs to be done before running the script.
        """

        owd = os.getcwd()
        # need to add absolute path here because otherwise bliding script doesnt find the file
        name_twopointfile = os.path.join(paths['runs_config'], 'measurement', '%s_twopointfile.fits'%name)
        print('We will bind this file:', name_twopointfile)
        filename_absolute_path = owd[:-3] + name_twopointfile[2:]
        os.chdir(path_2pt_pipeline)
        os.system('python pipeline/blind_2pt_usingcosmosis.py -s I_desperately_need_coffee -i pipeline/blinding_params_%s.ini -b add -u %s'%(config['lens_v'], filename_absolute_path))
        os.chdir(owd)
        if os.path.exists('%s_BLINDED.fits'%name_twopointfile[:-5]):
            os.system('rm %s'%name_twopointfile)

blind_measurements('gt')
blind_measurements('gt_boosted')
