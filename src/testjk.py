import numpy as np
import pyfits as pf
from jk import jk, plot_jk
from info import paths,config

"""
data = pf.getdata('../cats/y3/redmagic/v1/y3_gold_1.0_wide_cmv02a_cmcm2_run_redmapper_v6.4.20_redmagic_highdens_0.5-10_randoms.fit')
rand = np.random.random(len(data['ra']))
maskr = rand<0.001
ra = data['ra'][maskr]
dec = data['dec'][maskr]
jk = jk(paths,config,ra,dec)
"""

def get_lens(lens):
    """
    Given a lens sample, returns ra, dec, jk and weight, in case it exists.
    """

    ra_l = lens['ra']
    dec_l = lens['dec']
    jk_l = lens['jk']
    try:
        w_l = lens['w']
        print 'Weights found in lens catalog.'
    except:
        print 'There are no identified weights for the lenses.'
        w_l = np.ones(len(ra_l))

    return ra_l, dec_l, jk_l, w_l


lens_all = pf.getdata(paths['y1'] + 'lens.fits')
ra, dec, jk, w = get_lens(lens_all)
wjk = np.loadtxt('/Volumes/Data/y1_shear_tests/scripts/wrong_jks').astype('int')
print wjk
mask = np.in1d(jk,wjk)
print mask.sum()
plot_jk(paths, config, ra[mask], dec[mask], jk[mask])

