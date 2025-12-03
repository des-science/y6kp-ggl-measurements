import numpy as np
import h5py as h5
import joblib
import pickle
from astropy.io import fits
import healsparse as hsp

def load_lens_Y6_maglim(file_name, weights_file_name, zbin):
    """Load ra, dec, weight of Y6 MagLim++ zbin redshift bin from file"""
   
    with h5.File(file_name, 'r') as f:
        
        if zbin != 'all':  
            ra = np.array(f[f'desy6kp/maglim/tomo_bin_{zbin}/RA'])
            dec = np.array(f[f'desy6kp/maglim/tomo_bin_{zbin}/DEC'])
            w = np.array(f[f'desy6kp/maglim/tomo_bin_{zbin}/weight'])
        else:
            ra = []
            dec = []
            for zi in range(6):
                ra_i = np.array(f[f'desy6kp/maglim/tomo_bin_{zi}/RA'])
                dec_i = np.array(f[f'desy6kp/maglim/tomo_bin_{zi}/DEC'])
                w_i = np.array(f[f'desy6kp/maglim/tomo_bin_{zi}/weight'])
                ra.append(ra_i)
                dec.append(dec_i)
                w.append(w_i)
            ra = np.array(np.concatenate(ra)) 
            dec = np.array(np.concatenate(dec))
            w = np.array(np.concatenate(w))
    
            
    #if zbin != 'all':   
    #    w = joblib.load(weights_file_name)[zbin]
    #else:
    #    w = []
    #    for zi in range(6):
    #        w_i = joblib.load(weights_file_name)[zi]
    #        w.append(w_i)
    #    w = np.array(np.concatenate(w))
        
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
        if zbin != 'all':  
            ra = np.array(f[f'desy6kp/ran_maglim/tomo_bin_{zbin}/ra'])
            dec = np.array(f[f'desy6kp/ran_maglim/tomo_bin_{zbin}/dec'])
        else:
            ra = []
            dec = []
            for zbin in range(6):
                ra_i = np.array(f[f'desy6kp/ran/tomo_bin_{zbin}/ra'])
                dec_i = np.array(f[f'desy6kp/ran/tomo_bin_{zbin}/dec'])
                ra.append(ra_i)
                dec.append(dec_i)
            ra = np.array(np.concatenate(ra))
            dec = np.array(np.concatenate(dec))    

    return ra, dec
    

def load_randoms_Y3(file_name, zbin_lims):
    """Load ra, dec of Y6 Randoms in [zbin_lims[0], zbin_lims[1]] redshift bin from file"""
    
    with fits.open(file_name) as f:
        ra = f[1].data['RA']
        dec = f[1].data['DEC']
        z = f[1].data['Z']

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
    #with open(binning_file_name, 'rb') as file:
    #    tomo_bin = pickle.load(file)
    #    if zbin != 'all':
    #        mask_bin = np.where(tomo_bin==zbin)[0]
    #    else:
    #        mask_bin = np.where(tomo_bin!=-100)[0]

    with fits.open(file_name) as f:
        mask_bin = f[1].data['tomo']==zbin

        ra =  np.array(f[1].data['ra'])[mask_bin]
        dec = np.array(f[1].data['dec'])[mask_bin]
            
    mask = hsp.HealSparseMap.read(mask_file_name)
    mask_bool = mask.get_values_pos(ra, dec)
        
    ### Temporary: new BFD catalog - no need for logPQR()
    with fits.open(file_name) as f:
        # P = np.array(f[1].data['P'])[mask_bin][mask_bool]
        # Q0 = np.array(f[1].data['Q0'])[mask_bin][mask_bool]
        # Q1 = np.array(f[1].data['Q1'])[mask_bin][mask_bool]
        # R00 = np.array(f[1].data['R00'])[mask_bin][mask_bool]
        # R01 = np.array(f[1].data['R01'])[mask_bin][mask_bool]
        # R11 = np.array(f[1].data['R11'])[mask_bin][mask_bool]
        logpqr = np.array(f[1].data['pqr'])[mask_bin][mask_bool]
            
    ra =  ra[mask_bool]
    dec = dec[mask_bool]
            
    #pqr = ((np.vstack([P, Q0, Q1, R00, R01, R11])).T).astype(np.float64)
    #logpqr = logPQR(pqr)

    P = logpqr[:,0] 
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
        
        # Weighted averages of each of the terms to compute the response
        e1p = np.average(np.array(f[f'desy6kp/mdet/1p/tomo_bin_{zbin}/gauss_g_1']), 
                         weights = np.array(f[f'desy6kp/mdet/1p/tomo_bin_{zbin}/w'])) 
        e1m = np.average(np.array(f[f'desy6kp/mdet/1m/tomo_bin_{zbin}/gauss_g_1']), 
                         weights = np.array(f[f'desy6kp/mdet/1m/tomo_bin_{zbin}/w'])) 
        e2p = np.average(np.array(f[f'desy6kp/mdet/2p/tomo_bin_{zbin}/gauss_g_2']), 
                         weights = np.array(f[f'desy6kp/mdet/2p/tomo_bin_{zbin}/w'])) 
        e2m = np.average(np.array(f[f'desy6kp/mdet/2m/tomo_bin_{zbin}/gauss_g_2']), 
                         weights = np.array(f[f'desy6kp/mdet/2m/tomo_bin_{zbin}/w'])) 
    
    # Compute weighted average response
    deltag = 0.01 # shear applied to each image
    R11 = (e1p-e1m)/(2*deltag)
    R22 = (e2p-e2m)/(2*deltag)
    
    # Combine the components
    R = (R11+R22)/2.

    return ra, dec, e1, e2, R, w


def load_source_metacal(file_name, zbin):

    with h5.File(file_name, 'r') as f:
        
        select_metacal = np.array(f['index/select'])
        select_zbin = np.where(np.array(f['catalog/sompz/unsheared/bhat'])[select_metacal]==zbin)[0]

        ra = np.array(f[f'catalog/metacal/unsheared/ra'])[select_metacal][select_zbin]
        dec = np.array(f[f'catalog/metacal/unsheared/dec'])[select_metacal][select_zbin]
        e1 = np.array(f[f'catalog/metacal/unsheared/e_1'])[select_metacal][select_zbin]
        e2 = np.array(f[f'catalog/metacal/unsheared/e_2'])[select_metacal][select_zbin]
        w = np.array(f[f'catalog/metacal/unsheared/weight'])[select_metacal][select_zbin]
        
        select_metacal = np.array(f['index/select_1p'])
        select_zbin = np.where(np.array(f['catalog/sompz/sheared_1p/bhat'])[select_metacal]==zbin)[0]
        e1p = np.average(np.array(f[f'catalog/metacal/unsheared/e_1'])[select_metacal][select_zbin],
                         weights = np.array(f[f'catalog/metacal/sheared_1p/weight'])[select_metacal][select_zbin])
        
        select_metacal = np.array(f['index/select_1m'])
        select_zbin = np.where(np.array(f['catalog/sompz/sheared_1m/bhat'])[select_metacal]==zbin)[0]
        e1m = np.average(np.array(f[f'catalog/metacal/unsheared/e_1'])[select_metacal][select_zbin],
                         weights = np.array(f[f'catalog/metacal/sheared_1m/weight'])[select_metacal][select_zbin])
        
        select_metacal = np.array(f['index/select_2p'])
        select_zbin = np.where(np.array(f['catalog/sompz/sheared_2p/bhat'])[select_metacal]==zbin)[0]
        e2p = np.average(np.array(f[f'catalog/metacal/unsheared/e_2'])[select_metacal][select_zbin],
                         weights = np.array(f[f'catalog/metacal/sheared_2p/weight'])[select_metacal][select_zbin])
        
        select_metacal = np.array(f['index/select_2m'])
        select_zbin = np.where(np.array(f['catalog/sompz/sheared_2m/bhat'])[select_metacal]==zbin)[0]
        e2m = np.average(np.array(f[f'catalog/metacal/unsheared/e_2'])[select_metacal][select_zbin],
                         weights = np.array(f[f'catalog/metacal/sheared_2m/weight'])[select_metacal][select_zbin])
        
        deltag = 0.01 # shear applied to each image
        R11s = (e1p-e1m)/(2*deltag)
        R22s = (e2p-e2m)/(2*deltag)
        
        select_metacal = np.array(f['index/select'])
        select_zbin = np.where(np.array(f['catalog/sompz/unsheared/bhat'])[select_metacal]==zbin)[0]
        R11gs = np.array(f[f'catalog/metacal/unsheared/R11'])[select_metacal][select_zbin]
        R22gs = np.array(f[f'catalog/metacal/unsheared/R22'])[select_metacal][select_zbin]
        
        R11 = R11gs + R11s
        R22 = R22gs + R22s
        R11 = np.average(R11, weights=np.array(f[f'catalog/metacal/unsheared/weight'])[select_metacal][select_zbin])
        R22 = np.average(R22, weights=np.array(f[f'catalog/metacal/unsheared/weight'])[select_metacal][select_zbin])
        R = (R11 + R22)/2.
        
    return ra, dec, e1, e2, R, w

    