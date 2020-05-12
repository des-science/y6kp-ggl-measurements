import os

sample = 'redmagic'
#sample = 'maglim'

#scales = 'KP'
scales = 'HOD'

path_2pt_pipeline = '../../2pt_pipeline/'


if sample == 'redmagic':
        if scales == 'KP':
                file_to_blind = '/global/homes/j/jprat/y3-3x2pt-methods/cosmosis/data_vectors/v0.40_fiducial.fits'
                ini_ending = sample
        if scales == 'HOD':
                file_to_blind = '/global/homes/j/jprat/y3-3x2pt-methods/cosmosis/data_vectors/sim_HOD.fits'
                ini_ending = sample + '_' + scales

if sample == 'maglim':
        file_to_blind = '/global/homes/j/jprat/y3-3x2pt-methods/cosmosis/data_vectors/fiducial_maglim_cov_sourcesv040.fits'
        ini_ending = sample

print(file_to_blind)

def blind_measurements():
        """
        Run blinding script from 2pt pipeline repository. 
        In order to run this, cosmosis environment needs to be sourced.
        It needs to be done before running the script.
        """

        owd = os.getcwd()
        os.chdir(path_2pt_pipeline)
        print('We are using this params file: blinding_params_%s.ini'%ini_ending)
        os.system("python pipeline/blind_2pt_usingcosmosis.py -i pipeline/blinding_params_%s.ini -b add -u %s"%(ini_ending,file_to_blind))
        os.chdir(owd)

blind_measurements()
