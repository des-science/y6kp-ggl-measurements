import os

path_2pt_pipeline = '../../2pt_pipeline/'
file_to_blind = '/global/homes/j/jprat/y3-3x2pt-methods/cosmosis/data_vectors/v0.40_fiducial.fits'

def blind_measurements():
        """
        Run blinding script from 2pt pipeline repository. 
        In order to run this, cosmosis environment needs to be sourced.
        It needs to be done before running the script.
        """

        owd = os.getcwd()
        os.chdir(path_2pt_pipeline)
        os.system("python pipeline/blind_2pt_usingcosmosis.py -i pipeline/blinding_params_template_redmagic.ini -b add -u %s"%file_to_blind)
        os.chdir(owd)

blind_measurements()
