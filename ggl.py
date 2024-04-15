import numpy as np
import treecorr


def save_results(ggl_file_out, theta, gammat, gammax, gammat_shot_noise, xi_npairs, xi_weight,
                gammat_tot=None, gammax_tot=None, boost=None,
                gammat_rand=None, gammax_rand=None, xi_npairs_rand=None, xi_weight_rand=None):
    
    if ggl_file_out:
    
        if (boost is not None) & (gammat_rand is not None):
            np.savetxt(ggl_file_out, np.c_[theta, gammat_tot, gammat, gammax_tot, gammax, gammat_rand, gammax_rand, boost, gammat_shot_noise, xi_npairs, xi_npairs_rand, xi_weight, xi_weight_rand],
                        header=f'theta gammat_bf_rp gammat gammax_bf_rp gammax gammat_rand gammax_rand boost gammat_shot_noise xi_npairs xi_npairs_rand xi_weight xi_weight_rand\n')

        elif boost is not None:
            np.savetxt(ggl_file_out, np.c_[theta, gammat_tot, gammat, gammax_tot, gammax, boost, gammat_shot_noise, xi_npairs, xi_weight],
                        header=f'theta gammat_bf gammat gammax_bf gammax boost gammat_shot_noise xi_npairs xi_weight\n')

        elif gammat_rand is not None:
            np.savetxt(ggl_file_out, np.c_[theta, gammat_tot, gammat, gammax_tot, gammax, gammat_rand, gammax_rand, gammat_shot_noise, xi_npairs, xi_npairs_rand, xi_weight, xi_weight_rand],
                        header=f'theta gammat_rp gammat gammax_rp gammax gammat_rand gammax_rand gammat_shot_noise xi_npairs xi_npairs_rand xi_weight xi_weight_rand\n')

        else:
            np.savetxt(ggl_file_out, np.c_[theta, gammat, gammax, gammat_shot_noise, xi_npairs, xi_weight],
                        header=f'theta gammat gammax gammat_shot_noise xi_npairs xi_weight\n')

        print('GGL results saved')
        
    else:
        print('Not saving GGL results to file')
        
    return


def save_cov(cov_file_out, cov_gammat, cov_gammax, cov_boost=None, cov_gammat_rand=None):
    
    if cov_file_out:
    
        cov_gammat_file ='_gammat.'.join( f'{cov_file_out}'.rsplit('.',1))
        cov_gammax_file ='_gammax.'.join( f'{cov_file_out}'.rsplit('.',1))

        with open(cov_gammat_file,'wb') as f:
            for line in cov_gammat.T:
                np.savetxt(f, [line])      

        with open(cov_gammax_file, 'wb') as f:
            for line in cov_gammax.T:
                np.savetxt(f, [line])

        if cov_boost is not None:
            cov_boost_file ='_boost.'.join( f'{cov_file_out}'.rsplit('.',1))
            with open(cov_boost_file, 'wb') as f:
                for line in cov_boost.T:
                    np.savetxt(f, [line])

        if cov_gammat_rand is not None:       
            cov_gammat_rand_file ='_rand.'.join( f'{cov_file_out}'.rsplit('.',1))
            with open(cov_gammat_rand_file, 'wb') as f:
                for line in cov_gammat_rand.T:
                    np.savetxt(f, [line])

        print('Covariances saved')

    else:
        print('Not saving covariances to file')
        
    return


def boost_factor_calculate(sum_w_l, sum_w_r, w_LS, w_RS): 
    """Calculate boost factors as a function of theta
    - w_l: weights associated with lenses
    - w_r: weights associated with randoms
    - w_LS: weights associated with lens-source pairs
    - w_RS: weights associated with random-source pairs
    """
    
    boost = (sum_w_r/sum_w_l) * (w_LS/w_RS)
    return boost
        

def get_ggl(ra_l, dec_l, w_l, ra_s, dec_s, w_s, 
               e1, e2, Rg, ra_r, dec_r, units,
               theta_lims, nbins, sep_units, bin_slop, low_mem,
               use_randoms=False, use_boost=False, 
               compute_cov=False, npatch=None, ggl_file_out='', cov_file_out=''):
    """ Calculate GGL
    """
    
    print('Running GGL ...')
    # minimum and maximum angular separation
    theta_min, theta_max = theta_lims
    
    treecorr_config = {'nbins': nbins,
                        'min_sep': theta_min,
                        'max_sep': theta_max,
                        'sep_units': sep_units,
                        'bin_slop': bin_slop
                    }
    patch_centers = None
    
    if compute_cov:
        
        treecorr_config['var_method'] = 'jackknife'
        patch_centers='jk_centers'

        # generate catalog to save jk_centers 
        if ra_r is not None:
            cat_r = treecorr.Catalog(ra=ra_r, dec=dec_r, ra_units=units, dec_units=units, npatch=npatch)
            cat_r.write_patch_centers(patch_centers)
        else:
            ### seems ok from treecorr repo, or with sources ?
            cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l, ra_units=units, dec_units=units, npatch=npatch)
            cat_l.write_patch_centers(patch_centers)
        
    # generate lens and source catalogs to correlate and process them
    ng = treecorr.NGCorrelation(treecorr_config)

    cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l, ra_units=units, dec_units=units, w=w_l, patch_centers=patch_centers)
        
    cat_s = treecorr.Catalog(ra=ra_s, dec=dec_s, ra_units=units, dec_units=units,
                             g1=(e1-np.average(e1, weights=w_s)), g2=(e2-np.average(e2, weights=w_s)), 
                             w=w_s, patch_centers=patch_centers)
    
    ng.process(cat_l, cat_s, low_mem=low_mem)
    print('Done NG')

    # get theta, gammat
    theta = np.exp(ng.logr)
    gammat = ng.xi/Rg

    # get gammax
    gammax = ng.xi_im/Rg

    gammat_shot_noise = ng.varxi
    xi_npairs = ng.npairs
    xi_weight = ng.weight
    
    gammat_tot = None
    gammax_tot = None
    boost = None
    gammat_rand = None
    gammax_rand = None
    xi_npairs_rand = None
    xi_weight_rand = None
    
    if not (use_randoms or use_boost):
        save_results(ggl_file_out, theta, gammat, gammax, gammat_shot_noise, xi_npairs, xi_weight)

    elif (use_randoms or use_boost):
        
        if not compute_cov:
            # generate randoms catalogs to correlate and process them
            cat_r = treecorr.Catalog(ra=ra_r, dec=dec_r, ra_units=units, dec_units=units, patch_centers=patch_centers)

        rg = treecorr.NGCorrelation(treecorr_config)
        rg.process(cat_r, cat_s, low_mem=low_mem)
        print('Done RG')

        gammat_tot = np.copy(gammat)
        gammax_tot = np.copy(gammax)

        # get boost factors for gammat
        if use_boost:
            sum_w_l = np.sum(w_l)
            sum_w_r = len(ra_r) ### no weights for randoms ?
            boost = boost_factor_calculate(sum_w_l, sum_w_r, ng.weight, rg.weight)
            gammat_tot *= boost

        # random-point subtraction for gammat and gammax
        if use_randoms:
            gammat_rand = rg.xi/Rg
            gammax_rand = rg.xi_im/Rg
            xi_npairs_rand = rg.npairs
            xi_weight_rand = rg.weight
            gammat_tot -= gammat_rand
            gammax_tot -= gammax_rand

        assert not (gammat_tot == gammat).all()
        assert not (gammax_tot == gammax).all()
            
        save_results(ggl_file_out, theta, gammat, gammax, gammat_shot_noise, xi_npairs, xi_weight,
                    gammat_tot, gammax_tot, boost,
                    gammat_rand, gammax_rand, xi_npairs_rand, xi_weight_rand)

    cov_gammat = None
    cov_gammax = None
    cov_boost = None
    cov_gammat_rand = None
        
    # get gammat gammat covariance
    if compute_cov:
        
        # update correlations with responses to use in Jackknife mode
        ng.Rg = Rg*np.ones(len(theta))
            
        # generate fake treecorr correlation objects for lenses and randoms 
        # that hold the weights for the boost factor covariance calculations
        if use_boost:
            # initialize NN correlations for single point cross lenses and randoms
            nn_lp = treecorr.NNCorrelation(nbins=nbins, min_sep=1.e-3, max_sep=1.e5, 
                                            sep_units='arcmin', bin_slop=bin_slop, var_method='jackknife')
            nn_rp = treecorr.NNCorrelation(nbins=nbins, min_sep=1.e-3, max_sep=1.e5, 
                                            sep_units='arcmin', bin_slop=bin_slop, var_method='jackknife')
            # catalog containing single point
            cat_p = treecorr.Catalog(ra=np.array([np.mean(ra_l)]), dec=np.array([np.mean(dec_l)]), ra_units=units, dec_units=units)
            # process catalogs and get objects containing the weights split into jackknife patches
            nn_lp.process(cat_l, cat_p, low_mem=low_mem)
            nn_rp.process(cat_r, cat_p, low_mem=low_mem)
        
        # get gammat covariance
        if use_randoms:
            rg.Rg = Rg*np.ones(len(theta))

            if use_boost:
                func = lambda corrs: ( (corrs[0].weight/np.sum(corrs[2].weight)) / (corrs[1].weight/np.sum(corrs[3].weight)) * corrs[0].xi/corrs[0].Rg - corrs[1].xi/corrs[1].Rg )
                corrs = [ng,rg,nn_lp,nn_rp]
            else:
                func = lambda corrs: ( corrs[0].xi/corrs[0].Rg - corrs[1].xi/corrs[1].Rg )
                corrs = [ng,rg]
        else:
            if use_boost:
                func = lambda corrs: ( (corrs[0].weight/np.sum(corrs[2].weight)) / (corrs[1].weight/np.sum(corrs[3].weight)) 
                                        * corrs[0].xi/corrs[0].Rg )
                corrs = [ng,rg,nn_lp,nn_rp]
            else:
                func = lambda corrs: corrs[0].xi / corrs[0].Rg
                corrs = [ng]
        cov_gammat = treecorr.estimate_multi_cov(corrs, 'jackknife', func=func)

        # get gammax covariance
        if use_randoms:
            func = lambda corrs: corrs[0].xi_im/corrs[0].Rg - corrs[1].xi_im/corrs[1].Rg
            corrs = [ng,rg]
        else:
            func = lambda corrs: corrs[0].xi_im/corrs[0].Rg
            corrs = [ng]
        cov_gammax = treecorr.estimate_multi_cov(corrs, 'jackknife', func=func)

        # get boost factor covariance
        if use_boost:
            func = lambda corrs: (corrs[0].weight/np.sum(corrs[2].weight)) / (corrs[1].weight/np.sum(corrs[3].weight))
            corrs = [ng,rg,nn_lp,nn_rp]
            cov_boost = treecorr.estimate_multi_cov(corrs, 'jackknife', func=func)

        # get covariance of randoms points
        if use_randoms:
            func = lambda corrs: corrs[0].xi/corrs[0].Rg
            corrs = [rg]
            cov_gammat_rand = treecorr.estimate_multi_cov(corrs, 'jackknife', func=func)

        print('Done Covariance')

        save_cov(cov_file_out, cov_gammat, cov_gammax, cov_boost, cov_gammat_rand)
        
    assert not (gammat_tot == gammat).all()
    assert not (gammax_tot == gammax).all()

    return (theta, gammat, gammax, gammat_shot_noise, xi_npairs, xi_weight,
            gammat_tot, gammax_tot, boost,
            gammat_rand, gammax_rand, xi_npairs_rand, xi_weight_rand,
            cov_gammat, cov_gammax, cov_boost, cov_gammat_rand)


### maybe code better and have just one function get_ggl() for both mdet and bfd
def get_ggl_bfd(ra_l, dec_l, w_l, ra_s, dec_s, w_s, 
               Q0, Q1, R00, R01, R11, ra_r, dec_r, units,
               theta_lims, nbins, sep_units, bin_slop, low_mem,
               use_randoms=False, use_boost=False, 
               compute_cov=False, npatch=None, ggl_file_out='', cov_file_out=''):
    
    """Calculate GGL
    """  
    
    print('Running GGL ...')
    # minimum and maximum angular separation
    theta_min, theta_max = theta_lims
    
    treecorr_config = {'nbins': nbins,
                        'min_sep': theta_min,
                        'max_sep': theta_max,
                        'sep_units': sep_units,
                        'bin_slop': bin_slop
                    }
    patch_centers = None
    
    r = (R00 + R11)/2
    q = (R00 - R11)/2
        
    if compute_cov:
        
        treecorr_config['var_method'] = 'jackknife'
        patch_centers='jk_centers'

        # generate catalog to save jk_centers 
        if ra_r is not None:
            cat_r = treecorr.Catalog(ra=ra_r, dec=dec_r, ra_units=units, dec_units=units, npatch=npatch)
            cat_r.write_patch_centers(patch_centers)
        else:
            ### seems ok from treecorr repo, or with sources ?
            cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l, ra_units=units, dec_units=units, npatch=npatch)
            cat_l.write_patch_centers(patch_centers)
        
    # generate lens and source catalogs to correlate and process them
    cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l, ra_units=units, dec_units=units, w=w_l, patch_centers=patch_centers)
    
    cat_s = treecorr.Catalog(ra=ra_s, dec=dec_s, ra_units=units, dec_units=units, ### w=w_s ?
                            g1=Q0, g2=Q1, k=np.real(r), q1=np.real(q), q2=np.imag(q), patch_centers=patch_centers)
    
    ng = treecorr.NGCorrelation(treecorr_config)
    ng.process(cat_l, cat_s, low_mem=low_mem)
    nk = treecorr.NKCorrelation(treecorr_config)
    nk.process(cat_l, cat_s, low_mem=low_mem)
    nq = treecorr.NQCorrelation(treecorr_config)
    nq.process(cat_l, cat_s, low_mem=low_mem)
    print('Done NG')
    
    r = nk.xi 
    q = nq.xi + 1j * nq.xi_im
    Q = ng.xi + 1j * ng.xi_im
    g = (r * Q - q * np.conj(Q)) / (np.abs(r)**2 - np.abs(q)**2)
    
    # get theta, gammat
    theta = np.exp(ng.logr)
    gammat = np.real(g)

    # get gammax
    gammax = np.imag(g)

    gammat_shot_noise = ng.varxi ### not sure at all
    xi_npairs = ng.npairs ### maybe
    xi_weight = ng.weight ### not sure --> not sure about boost
    
    gammat_tot = None
    gammax_tot = None
    boost = None
    gammat_rand = None
    gammax_rand = None
    xi_npairs_rand = None
    xi_weight_rand = None

    if not (use_randoms or use_boost):
        save_results(ggl_file_out, theta, gammat, gammax, gammat_shot_noise, xi_npairs, xi_weight)

    elif (use_randoms or use_boost):
        
        if not compute_cov:
            # generate randoms catalogs to correlate and process them
            cat_r = treecorr.Catalog(ra=ra_r, dec=dec_r, ra_units=units, dec_units=units, patch_centers=patch_centers)
        
        rg = treecorr.NGCorrelation(treecorr_config)
        rg.process(cat_r, cat_s, low_mem=low_mem)
        rk = treecorr.NKCorrelation(treecorr_config)
        rk.process(cat_r, cat_s, low_mem=low_mem)
        rq = treecorr.NQCorrelation(treecorr_config)
        rq.process(cat_r, cat_s, low_mem=low_mem)
        print('Done RG')
        
        r_rg = rk.xi 
        q_rg = rq.xi + 1j * rq.xi_im
        Q_rg = rg.xi + 1j * rg.xi_im
        g_rg = (r_rg * Q_rg - q_rg * np.conj(Q_rg)) / (np.abs(r_rg)**2 - np.abs(q_rg)**2)

        gammat_tot = np.copy(gammat)
        gammax_tot = np.copy(gammax)
        
        # get boost factors for gammat
        if use_boost:
            sum_w_l = np.sum(w_l)
            sum_w_r = len(ra_r) ### no weights for randoms ?
            boost = boost_factor_calculate(sum_w_l, sum_w_r, ng.weight, rg.weight)
            gammat_tot *= boost
            
        # random-point subtraction for gammat and gammax
        if use_randoms:
            gammat_rand = np.real(g_rg)
            gammax_rand = np.imag(g_rg)
            xi_npairs_rand = rg.npairs ### maybe
            xi_weight_rand = rg.weight  ### not sure
            gammat_tot -= gammat_rand
            gammax_tot -= gammax_rand

        assert not (gammat_tot == gammat).all()
        assert not (gammax_tot == gammax).all()
            
        save_results(ggl_file_out, theta, gammat, gammax, gammat_shot_noise, xi_npairs, xi_weight,
                    gammat_tot, gammax_tot, boost,
                    gammat_rand, gammax_rand, xi_npairs_rand, xi_weight_rand)

    cov_gammat = None
    cov_gammax = None
    cov_boost = None
    cov_gammat_rand = None
    
    # get gammat gammat covariance
    ### try to find a better way to code this
    if compute_cov:
        
        # update correlations with responses to use in Jackknife mode
        ### ng.Rg = Rg*np.ones(len(theta)) # no Rg in BFD ?
            
        # generate fake treecorr correlation objects for lenses and randoms 
        # that hold the weights for the boost factor covariance calculations
        if use_boost:
            # initialize NN correlations for single point cross lenses and randoms
            nn_lp = treecorr.NNCorrelation(nbins=nbins, min_sep=1.e-3, max_sep=1.e5, 
                                            sep_units='arcmin', bin_slop=bin_slop, var_method='jackknife')
            nn_rp = treecorr.NNCorrelation(nbins=nbins, min_sep=1.e-3, max_sep=1.e5, 
                                            sep_units='arcmin', bin_slop=bin_slop, var_method='jackknife')
            # catalog containing single point
            cat_p = treecorr.Catalog(ra=np.array([np.mean(ra_l)]), dec=np.array([np.mean(dec_l)]), ra_units=units, dec_units=units)
            # process catalogs and get objects containing the weights split into jackknife patches
            nn_lp.process(cat_l, cat_p, low_mem=low_mem)
            nn_rp.process(cat_r, cat_p, low_mem=low_mem)
        
        # get gammat covariance
        ### are funcs ok ?
        if use_randoms:
            ### rg.Rg = Rg*np.ones(len(theta)) # no Rg in BFD ?

            if use_boost:
                func = lambda corrs: (corrs[0].weight/np.sum(corrs[6].weight)) / (corrs[3].weight/np.sum(corrs[7].weight)) * np.real((corrs[1].xi * corrs[0].xi + 1j * corrs[0].xi_im - corrs[2].xi + 1j * corrs[2].xi_im * np.conj(corrs[0].xi + 1j * corrs[0].xi_im)) / (np.abs(corrs[1].xi)**2 - np.abs(corrs[2].xi + 1j * corrs[2].xi_im)**2)) - np.real((corrs[4].xi * corrs[3].xi + 1j * corrs[3].xi_im - corrs[5].xi + 1j * corrs[5].xi_im * np.conj(corrs[3].xi + 1j * corrs[3].xi_im)) / (np.abs(corrs[4].xi)**2 - np.abs(corrs[5].xi + 1j * corrs[5].xi_im)**2))
                corrs = [ng, nk, nq, rg, rk, rq, nn_lp, nn_rp]
            else:
                func = lambda corrs: np.real((corrs[1].xi * corrs[0].xi + 1j * corrs[0].xi_im - corrs[2].xi + 1j * corrs[2].xi_im * np.conj(corrs[0].xi + 1j * corrs[0].xi_im)) / (np.abs(corrs[1].xi)**2 - np.abs(corrs[2].xi + 1j * corrs[2].xi_im)**2)) - np.real((corrs[4].xi * corrs[3].xi + 1j * corrs[3].xi_im - corrs[5].xi + 1j * corrs[5].xi_im * np.conj(corrs[3].xi + 1j * corrs[3].xi_im)) / (np.abs(corrs[4].xi)**2 - np.abs(corrs[5].xi + 1j * corrs[5].xi_im)**2))
                corrs = [ng, nk, nq, rg, rk, rq]
        else:
            if use_boost:
                func = lambda corrs: (corrs[0].weight/np.sum(corrs[6].weight)) / (corrs[3].weight/np.sum(corrs[7].weight)) * np.real((corrs[1].xi * corrs[0].xi + 1j * corrs[0].xi_im - corrs[2].xi + 1j * corrs[2].xi_im * np.conj(corrs[0].xi + 1j * corrs[0].xi_im)) / (np.abs(corrs[1].xi)**2 - np.abs(corrs[2].xi + 1j * corrs[2].xi_im)**2))
                corrs = [ng, nk, nq, rg, rk, rq, nn_lp, nn_rp]
            else:
                func = lambda corrs: np.real((corrs[1].xi * corrs[0].xi + 1j * corrs[0].xi_im - corrs[2].xi  + 1j * corrs[2].xi_im * np.conj(corrs[0].xi + 1j * corrs[0].xi_im)) / (np.abs(corrs[1].xi)**2 - np.abs(corrs[2].xi + 1j * corrs[2].xi_im)**2))
                corrs = [ng, nk, nq]
        cov_gammat = treecorr.estimate_multi_cov(corrs, 'jackknife', func=func)

        # get gammax covariance
        if use_randoms:
            func = lambda corrs: np.imag((corrs[1].xi * corrs[0].xi + 1j * corrs[0].xi_im - corrs[2].xi + 1j * corrs[2].xi_im * np.conj(corrs[0].xi + 1j * corrs[0].xi_im)) / (np.abs(corrs[1].xi)**2 - np.abs(corrs[2].xi + 1j * corrs[2].xi_im)**2)) - np.imag((corrs[4].xi * corrs[3].xi + 1j * corrs[3].xi_im - corrs[5].xi + 1j * corrs[5].xi_im * np.conj(corrs[3].xi + 1j * corrs[3].xi_im)) / (np.abs(corrs[4].xi)**2 - np.abs(corrs[5].xi + 1j * corrs[5].xi_im)**2))
            corrs = [ng, nk, nq, rg, rk, rq]
        else:
            func = lambda corrs: np.imag((corrs[1].xi * corrs[0].xi + 1j * corrs[0].xi_im - corrs[2].xi + 1j * corrs[2].xi_im * np.conj(corrs[0].xi + 1j * corrs[0].xi_im)) / (np.abs(corrs[1].xi)**2 - np.abs(corrs[2].xi + 1j * corrs[2].xi_im)**2))
            corrs = [ng, nk, nq]
        cov_gammax = treecorr.estimate_multi_cov(corrs, 'jackknife', func=func)
        
        # get boost factor covariance
        if use_boost:
            func = lambda corrs: (corrs[0].weight/np.sum(corrs[6].weight)) / (corrs[3].weight/np.sum(corrs[7].weight))
            corrs = [ng, nk, nq, rg, rk, rq, nn_lp, nn_rp]
            cov_boost = treecorr.estimate_multi_cov(corrs, 'jackknife', func=func)
            
        # get covariance of randoms points
        if use_randoms:
            func = lambda corrs: np.real((corrs[1].xi * corrs[0].xi + 1j * corrs[0].xi_im - corrs[2].xi + 1j * corrs[2].xi_im * np.conj(corrs[0].xi + 1j * corrs[0].xi_im)) / (np.abs(corrs[1].xi)**2 - np.abs(corrs[2].xi + 1j * corrs[2].xi_im)**2))
            corrs = [rg, rk, rq]
            cov_gammat_rand = treecorr.estimate_multi_cov(corrs, 'jackknife', func=func)

        print('Done Covariance')

        save_cov(cov_file_out, cov_gammat, cov_gammax, cov_boost, cov_gammat_rand)
        
    assert not (gammat_tot == gammat).all()
    assert not (gammax_tot == gammax).all()

    return (theta, gammat, gammax, gammat_shot_noise, xi_npairs, xi_weight,
            gammat_tot, gammax_tot, boost,
            gammat_rand, gammax_rand, xi_npairs_rand, xi_weight_rand,
            cov_gammat, cov_gammax, cov_boost, cov_gammat_rand)
