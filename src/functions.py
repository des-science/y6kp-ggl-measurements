from __future__ import division
import numpy as np
#import cosmology
#from iminuit import Minuit

def inv_sigma_crit_eff(zlbin, zsbin, nzl, nzs, omegam = 0.3):
    '''
    Returns the inverse sigma_crit effective given a lens and source N(z)'s
    '''
    dzl = zlbin[1]-zlbin[0] 
    dzs = zsbin[1]-zsbin[0]
    norml = np.sum(nzl*dzl)
    norms = np.sum(nzs*dzs)
    nzl = nzl/norml
    nzs = nzs/norms
    isc = 0.
    norm = 0.

    # Define meshgrid for redshifts and for Nzs
    X,Y = np.meshgrid(zlbin, zsbin)
    NZL, NZS = np.meshgrid(nzl, nzs)
    # Construct 2-D integrand
    c=cosmology.Cosmo(omegam)
    sci_flat = c.sigmacritinv(X,Y)
    sci_re = np.reshape(sci_flat, (len(zsbin),len(zlbin)), order='C')
    integrand = NZL*NZS*sci_re
    # Do a 1-D integral over every row
    I = np.zeros(len(zsbin))
    for i in range(len(zsbin)):
        I[i] = np.trapz(integrand[i,:], zlbin)

    # Then an integral over the result
    F = np.trapz(I, zsbin)

    return F

def covariance(gt_all, gt_mean):
    '''
    Returns the covariance from many jackknife realizations.
    - gt_all: 2-D array of measurements with len(gt_all) = njk (number of jk regions) 
    - gt_mean: 1-D array with len(gt_mean) =  number of angular bins 
    '''
    cov = []
    for i in range(len(gt_mean)):
        r = []
        for l in range(len(gt_mean)):
            s = sum([ (gt_all[k][i] - gt_mean[i]) * (gt_all[k][l] - gt_mean[l]) for k in range(len(gt_all)) ])
            r.append(s)
        cov.append(r)

    COV = np.array(cov)
    COV = float(len(gt_all)-1)/len(gt_all) * COV
    print 'len(gt_all) in covariance', len(gt_all)
    return COV


def minimize_chi2(COV, data):
    INVCOV = np.linalg.inv(COV)
    # matrix cast operator
    data_mat = np.mat(data)
    # minimize the chi-square
    def f(data_fit):
        return (data_mat - data_fit)*INVCOV*(data_mat.T - data_fit)
    m = Minuit(f, print_level=0, errordef=1, pedantic = False)
    m.migrad()

    fit, err_fit = m.values['data_fit'], m.errors['data_fit']
    chi2_fit = f(fit)

    #print 'fit, err_fit = %0.3e +- %0.3e'%(fit, err_fit)
    #print 'chi2_fit/ndf: %0.2f/%d'%(chi2_fit[0], (len(data)-1))
    return fit, err_fit, chi2_fit[0], len(data)-1

def minimize_chi2_fit_amplitude(COV, data, sims):
    INVCOV = np.linalg.inv(COV)
    # matrix cast operator
    data_mat = np.mat(data)
    sims_mat = np.mat(sims)
    # minimize the chi-square
    def f(A):
        # No need to include A**(-2) because it's a constant and we are minimizing the chi2
        return (data_mat - A*sims_mat)*INVCOV*(data_mat.T - A*sims_mat.T) 
    def chi2(A):
        # To obtain the actual value for the chi2 yes
        return (data_mat - A*sims_mat)*INVCOV*(A**(-2))*(data_mat.T - A*sims_mat.T)
    m = Minuit(f, print_level=0, errordef=1, pedantic = False)
    m.migrad()

    fit, err_fit = m.values['A'], m.errors['A']
    chi2_fit = chi2(fit)

    return fit, chi2_fit[0], len(data)-1


def fit_constant(COV, d):
    # COV: covariance matrix
    # d: data vector to fit

    # Fit a constant
    mcov_inv = np.linalg.inv(COV)
    const = np.linspace(-0.001, 0.001, 1000001)
    chi2 = np.zeros(len(const))
    for i in range(len(const)):
        dmm = d - const[i]
        dmm = np.mat(dmm)
        #temp = dmm * mcov_inv
        #chi2[i] = np.sum(dmm * np.array(temp))
        chi2[i] = dmm * mcov_inv * dmm.T
        
    #print 'chi2:', chi2
    like = np.exp(-chi2 / 2.)
    #print 'like:', like
    tot = np.sum(like)
    #print 'tot:', tot
    prob = like / tot
    #print 'prob:', prob

    # Calculate cumulative distribution to read off confidence levels
    cum_arr = np.zeros(len(prob))
    cum = 0.
    for i in range(len(prob)):
        cum = cum + prob[i]
        cum_arr[i] = cum
    #print 'cum:', cum_arr

    # Find peak and mean of distribution
    ind = np.argmax(like)
    peak = const[ind]
    chi2_fit_peak = chi2[ind]
    ind = np.where(cum_arr < 0.5)[0]
    middle = const[ind[-1]]
    chi2_fit_middle = chi2[ind[-1]]
    
    # Find left and right edge of 1-sigma region
    ind = np.where(cum_arr < 0.15865)[0]
    left = const[ind[-1]]
    ind = np.where(cum_arr < 0.84135)[0]
    right = const[ind[-1]]

    std_dev = (right - left) / 2.

    #print '\n%s' % (zlmax[m])
    #print '\nleft, peak, right = %.2f, %.2f, %.2f' % (left, peak, right), '\n'
    #print 'chi2 null = ', chi2_null, ' / ', nR - nlow - nhi
    #print 'Number of constants: %d'%len(const)
    #print 'left, middle, right = %.3e, %.3e, %.3e' % (left, middle, right)
    #print 'left, peak, right = %.3e, %.3e, %.3e' % (left, peak, right)
    #print 'chi2 fit middle/ndf = %0.2f/%d'%(chi2_fit_middle, len(d)-1)
    #print 'chi2 fit peak/ndf = %0.2f/%d'%(chi2_fit_peak, len(d)-1)

    #if peak < 0:
    #    print 'Peak + right: ', peak + abs(right)
        
    #if peak > 0: 
    #    print 'Peak - left:', peak - abs(left)
    return  peak, left, right

def compute_sigmas(q1, err1, q2, err2):
    # q1, err1: quantity and corresponding error
    print 'q1,err1,q2,err2:', q1,err1,q2,err2
    return np.abs(q1-q2)/np.sqrt(err1**2 + err2**2)


