import numpy as np
import h5py as h5
import joblib
import pickle
from astropy.io import fits
import healsparse as hsp

def load_lens_Y6_maglim(file_name, weights_file_name, zbin):
    """Load ra, dec, weight of Y6 MagLim++ zbin redshift bin from file"""
                   
    with h5.File(file_name, 'r') as f:
        ra = np.array(f[f'desy6kp/maglim/tomo_bin_{zbin}/RA'])
        dec = np.array(f[f'desy6kp/maglim/tomo_bin_{zbin}/DEC'])
            
    w = joblib.load(weights_file_name)[zbin]
        
    return ra, dec, w
    

def load_lens_Y3_maglim(file_name, zbin_lims):
    """Load ra, dec, weight of Y3 MagLim redshift bin in [zbin_lims[0], zbin_lims[1]] from file"""
    
    with h5.File(file_name, 'r') as f:
        select_maglim = np.array(f['index/maglim/select'])
        ra = np.array(f['catalog/gold/ra'])[select_maglim]
        dec = np.array(f['catalog/gold/dec'])[select_maglim]
        z = np.array(f['catalog/dnf/unsheared/zmean_sof'])[select_maglim]
        
    ### Y3 maglim weights ?
    w = np.ones(len(ra))

    zbin = np.where((z>zbin_lims[0]) * (z<zbin_lims[1]))[0]
        
    return ra[zbin], dec[zbin], w[zbin]


def load_randoms_Y6(file_name, zbin):
    """Load ra, dec of Y6 Randoms zbin redshift bin from file"""
                                                
    with h5.File(file_name, 'r') as f:
        ra = np.array(f[f'desy6kp/ran/tomo_bin_{zbin}/RA'])
        dec = np.array(f[f'desy6kp/ran/tomo_bin_{zbin}/DEC'])

    ### random weights ?

    return ra, dec
    

def load_randoms_Y3(file_name, zbin_lims):
    """Load ra, dec of Y6 Randoms in [zbin_lims[0], zbin_lims[1]] redshift bin from file"""
    
    with fits.open(file_name) as f:
        ra = f[1].data['RA']
        dec = f[1].data['DEC']
        z = f[1].data['Z']

    ### random weights ?

    zbin = np.where((z>zbin_lims[0]) * (z<zbin_lims[1]))[0]
        
    return ra[zbin], dec[zbin]
    
    
def load_source_bfd(file_name, binning_file_name, mask_file_name, zbin):
    """Load Y6 BFD zbin redshift bin from file"""
     
        ### Temporary: : the bins assignment in the mastercat is wrong and e1,e2 in the bfd cat are wrong
        #with h5.File(file_name, 'r') as f:
        #    ra = np.array(f[f'desy6kp/bfd/tomo_bin_{zbin}/ra'])
        #    dec = np.array(f[f'desy6kp/bfd/tomo_bin_{zbin}/dec'])
        #    e1 = np.array(f[f'desy6kp/bfd/tomo_bin_{zbin}/e1'])
        #    e2 = np.array(f[f'desy6kp/bfd/tomo_bin_{zbin}/e2'])
        #    P = np.array(f[f'desy6kp/bfd/tomo_bin_{zbin}/P'])
        #    Q0 = np.array(f[f'desy6kp/bfd/tomo_bin_{zbin}/Q0'])
        #    Q1 = np.array(f[f'desy6kp/bfd/tomo_bin_{zbin}/Q1'])
        #    R00 = np.array(f[f'desy6kp/bfd/tomo_bin_{zbin}/R00'])
        #    R01 = np.array(f[f'desy6kp/bfd/tomo_bin_{zbin}/R01'])
        #    R11 = np.array(f[f'desy6kp/bfd/tomo_bin_{zbin}/R11'])
        
    with open(binning_file_name, 'rb') as file:
        tomo_bin = pickle.load(file)
    mask_bin = np.where(tomo_bin==zbin)[0]
            
    with fits.open(file_name) as f:
        ra =  np.array(f[1].data['ra'])[mask_bin]
        dec = np.array(f[1].data['dec'])[mask_bin]
            
    mask = hsp.HealSparseMap.read(mask_file_name)
    mask_bool = mask.get_values_pos(ra, dec)
        
    with fits.open(file_name) as f:
        P = np.array(f[1].data['P'])[mask_bin][mask_bool]
        Q0 = np.array(f[1].data['Q0'])[mask_bin][mask_bool]
        Q1 = np.array(f[1].data['Q1'])[mask_bin][mask_bool]
        R00 = np.array(f[1].data['R00'])[mask_bin][mask_bool]
        R01 = np.array(f[1].data['R01'])[mask_bin][mask_bool]
        R11 = np.array(f[1].data['R11'])[mask_bin][mask_bool]
            
    ra =  ra[mask_bool]
    dec = dec[mask_bool]
            
    pqr = ((np.vstack([P, Q0, Q1, R00, R01, R11])).T).astype(np.float64)
    logpqr = logPQR(pqr)

    Q0 = logpqr[:,1]
    Q1 = logpqr[:,2]
    R00 = logpqr[:,3]
    R01 = logpqr[:,4]
    R11 = logpqr[:,5]
            
    e1, e2 = approx_e_bfd(np.copy(logpqr))
        
    ### no weights in BFD ?
    w = np.ones(len(ra))
        
    return ra, dec, e1, e2, P, Q0, Q1, R00, R01, R11, w
    
    
def logPQR(pqr):
    """Convert an Nx6 array representing quadratic Taylor expansion of p w.r.t. shear into Taylor expansion of log(P).
        Any input rows with p<=0 are put as zeros in output
    """
    pqr_out = np.zeros_like(pqr)
    p = pqr[:,0]
    good_p = p>0.
    q = pqr[:,1:3] / p[:,np.newaxis]
    r = pqr[:,3:] / p[:,np.newaxis]
    pqr_out[:,0] = np.log(p)
    pqr_out[:,1:3] = q
    r[:,0] -= q[:,0]*q[:,0]
    r[:,1] -= q[:,1]*q[:,0]
    r[:,2] -= q[:,1]*q[:,1]
    pqr_out[:,3:] = r
    pqr_out = np.where(good_p[:,np.newaxis], pqr_out, 0.)
    return pqr_out
    
    
def approx_e_bfd(logpqr):
    """Compute BFD e1, e2 using the average response
    """
    q = logpqr[:,1:3]
    r = np.sum(logpqr[:,3:], axis=0, dtype=float)
    R = np.array([[r[0], r[1]], [r[1], r[2]]])
    invR = np.linalg.inv(R)
    q *= (len(q[:,0]))
    e = np.matmul(invR, q.T)  
    return e[0], e[1]
    

def load_source_metadetect(file_name, zbin):
    """Load Y6 Metadetect zbin redshift bin from file"""
            
    with h5.File(file_name, 'r') as f:
        ra = np.array(f[f'desy6kp/mdet/noshear/tomo_bin_{zbin}/ra'])
        dec = np.array(f[f'desy6kp/mdet/noshear/tomo_bin_{zbin}/dec'])
        e1 = np.array(f[f'desy6kp/mdet/noshear/tomo_bin_{zbin}/gauss_g_1'])
        e2 = np.array(f[f'desy6kp/mdet/noshear/tomo_bin_{zbin}/gauss_g_2'])
        w = np.array(f[f'desy6kp/mdet/noshear/tomo_bin_{zbin}/w'])
       
    R = compute_response(file_name, zbin)

    return ra, dec, e1, e2, R, w
    
def compute_response(file_name, zbin, deltag = 0.01):
    """Compute Metadetection weighted average response in the zbin redshift bin
    - deltag: shear applied to each image
    """
    
    with h5.File(file_name, 'r') as f:
    # Weighted averages of each of the terms
        e1p = np.average(f[f'desy6kp/mdet/1p/tomo_bin_{zbin}/gauss_g_1'], weights = f[f'desy6kp/mdet/1p/tomo_bin_{zbin}/w']) 
        e1m = np.average(f[f'desy6kp/mdet/1m/tomo_bin_{zbin}/gauss_g_1'], weights = f[f'desy6kp/mdet/1m/tomo_bin_{zbin}/w']) 
        e2p = np.average(f[f'desy6kp/mdet/2p/tomo_bin_{zbin}/gauss_g_2'], weights = f[f'desy6kp/mdet/2p/tomo_bin_{zbin}/w']) 
        e2m = np.average(f[f'desy6kp/mdet/2m/tomo_bin_{zbin}/gauss_g_2'], weights = f[f'desy6kp/mdet/2m/tomo_bin_{zbin}/w']) 
    
    # Compute the average response for each component
    R11 = (e1p-e1m)/(2*deltag)
    R22 = (e2p-e2m)/(2*deltag)
    
    # Combine the components
    R = (R11+R22)/2.

    return R
    
    """"
    TO DO
    def load_source_metacal(self, data_file, zbin=None):
        
        # maybe change it so that it reads the mastercat without calling destest
        
        # read yaml file that defines all the catalog selections used
        params = yaml.load(open(data_file), Loader=yaml.FullLoader)
        params['param_file'] = data_file.split('/')[-1]

        source_selector, source_calibrator = catalog_utils.load_catalog(params, 'mcal', 'mcal', 
                                                                        params['source_group'], params['source_table'], 
                                                                        params['source_path'], 
                                                                        return_calibrator=destest.MetaCalib)

        gold_selector = catalog_utils.load_catalog(params, 'gold', 'mcal', params['gold_group'], params['gold_table'],
                                                    params['gold_path'], inherit=source_selector)

        pz_selector = catalog_utils.load_catalog(params, 'pz', 'mcal', params['pz_group'], 
                                                    params['pz_table'], params['pz_path'], inherit=source_selector)

        # dictionary with the unsheared version of each quantity with the 
        # selections from: unsheared, 1p, 1m, 2p, 2m
        source_5sels = {}
        source_5sels['unsheared'] = {}
        source_5sels['unsheared']['ra'] = [gold_selector.get_col('ra', uncut=True)[0][gold_selector.get_mask()[i]] for i in range(5)]
        source_5sels['unsheared']['dec'] = [gold_selector.get_col('dec', uncut=True)[0][gold_selector.get_mask()[i]] for i in range(5)]
        source_5sels['unsheared']['e1'] = [source_selector.get_col('e_1', uncut=True)[0][source_selector.get_mask()[i]] for i in range(5)]
        source_5sels['unsheared']['e2'] = [source_selector.get_col('e_2', uncut=True)[0][source_selector.get_mask()[i]] for i in range(5)]

        # dictionary with the 5 selections (1p, 1m, 2p, 2m), in which the 
        # quantities are obtained from sheared images or using fluxes measured
        # on sheared images
        source_5sels['sheared'] = {}
        source_5sels['sheared']['e1'] = source_selector.get_col('e_1')
        source_5sels['sheared']['e2'] = source_selector.get_col('e_2')
        # formerly z_g
        source_5sels['sheared']['som_bin'] = pz_selector.get_col('bhat')

        # dictionary with the unsheared version and selection only
        source = {}
        source['ra'] = source_5sels['unsheared']['ra'][0]
        source['dec'] = source_5sels['unsheared']['dec'][0]
        source['e1'] = source_5sels['unsheared']['e1'][0]
        source['e2'] = source_5sels['unsheared']['e2'][0]
        source['som_bin'] = source_5sels['sheared']['som_bin'][0]

        # masks from redshifts
        photoz_masks = [
                        (source_5sels['sheared']['som_bin'][i] >= zs_bin) & \
                        (source_5sels['sheared']['som_bin'][i] < zs_bin+1) for i in range(5)
                       ]
        maskz = photoz_masks[0]

        ra_s = source['ra'][maskz]
        dec_s = source['dec'][maskz]
        e1_s = source['e1'][maskz]
        e2_s = source['e2'][maskz]

        # calibration given the photoz bin
        R1,_,w_g = source_calibrator.calibrate('e_1', mask=photoz_masks)
        R2,_,w_g = source_calibrator.calibrate('e_2', mask=photoz_masks)
        R_g = 0.5*(R1+R2)

        return ra_s, dec_s, e1_s, e2_s, R_g, w_g
"""
