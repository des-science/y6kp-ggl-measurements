import os
import matplotlib.pyplot as plt  
from info import blind, paths, config, zbins, plotting, source_nofz_pars, sysmaps, mode, filename_mastercat


path_2pt_pipeline = '../../2pt_pipeline/'

name_twopointfile = os.path.join(paths['runs_config'], 'measurement', 'gammat_twopointfile.fits')

def blind_measurements():
        """
        Run blinding script from 2pt pipeline repository. 
        In order to run this, cosmosis environment needs to be sourced.
        It needs to be done before running the script.
        """

        owd = os.getcwd()
        # need to add absolute path here because otherwise bliding script doesnt find the file
        filename_absolute_path = owd[:-3] + name_twopointfile[2:]
        os.chdir(path_2pt_pipeline)
        os.system('python pipeline/blind_2pt_usingcosmosis.py -s Y3_blinded -i pipeline/blinding_params_template.ini -b add -u %s'%filename_absolute_path)
        os.chdir(owd)
        if os.path.exists('%s_BLINDED.fits'%name_twopointfile[:-5]):
            os.system('rm %s'%name_twopointfile)

blind_measurements()
