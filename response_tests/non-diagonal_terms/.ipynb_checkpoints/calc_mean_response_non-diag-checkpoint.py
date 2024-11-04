import numpy as np
import h5py as h5

data_source = '/global/cfs/cdirs/des/giannini/ggl/data/2024-08-26/metadetect_v6_UNBLINDED_2024-08-26.hdf5'
mdet = h5.File(data_source, 'r')
path = '/global/cfs/projectdirs/des/awhyley/ggl/response_tests/non-diagonal_terms'

g_all = {}
w_all = {}

#Individual z bins

for zs_bin in range(4):

    print (zs_bin)

    g1p = np.array(mdet['1p']['tomo_bin_{}'.format(zs_bin)]['gauss_g_2'])
    g1m = np.array(mdet['1m']['tomo_bin_{}'.format(zs_bin)]['gauss_g_2'])
    w1p = np.array(mdet['1p']['tomo_bin_{}'.format(zs_bin)]['w'])
    w1m = np.array(mdet['1m']['tomo_bin_{}'.format(zs_bin)]['w'])
    R12 = (np.average(g1p, weights=w1p)-np.average(g1m, weights=w1m))/2/0.01

    g2p = np.array(mdet['2p']['tomo_bin_{}'.format(zs_bin)]['gauss_g_1'])
    g2m = np.array(mdet['2m']['tomo_bin_{}'.format(zs_bin)]['gauss_g_1'])
    w2p = np.array(mdet['2p']['tomo_bin_{}'.format(zs_bin)]['w'])
    w2m = np.array(mdet['2m']['tomo_bin_{}'.format(zs_bin)]['w'])
    R21 = (np.average(g2p, weights=w2p)-np.average(g2m, weights=w2m))/2/0.01
    
    response_components = np.array([R12, R21])
    np.savetxt(path+'/Response_non-diag_bin{}.txt'.format(zs_bin), response_components)
    
    #Now save the arrays to g_all so means from full catalogue can be computed
    
    g_all['g1p_bin{}'.format(zs_bin)] = g1p
    g_all['g1m_bin{}'.format(zs_bin)] = g1m
    g_all['g2p_bin{}'.format(zs_bin)] = g2p
    g_all['g2m_bin{}'.format(zs_bin)] = g2m
    
    w_all['w1p_bin{}'.format(zs_bin)] = w1p
    w_all['w1m_bin{}'.format(zs_bin)] = w1m
    w_all['w2p_bin{}'.format(zs_bin)] = w2p
    w_all['w2m_bin{}'.format(zs_bin)] = w2m
    
#Full catalogue

print("averaging ellipticities over all bins")

g1p_all = np.concatenate((g_all['g1p_bin0'], g_all['g1p_bin1'], g_all['g1p_bin2'], g_all['g1p_bin3']))
g1m_all = np.concatenate((g_all['g1m_bin0'], g_all['g1m_bin1'], g_all['g1m_bin2'], g_all['g1m_bin3']))
g2p_all = np.concatenate((g_all['g2p_bin0'], g_all['g2p_bin1'], g_all['g2p_bin2'], g_all['g2p_bin3']))
g2m_all = np.concatenate((g_all['g2m_bin0'], g_all['g2m_bin1'], g_all['g2m_bin2'], g_all['g2m_bin3']))

w1p_all = np.concatenate((w_all['w1p_bin0'], w_all['w1p_bin1'], w_all['w1p_bin2'], w_all['w1p_bin3']))
w1m_all = np.concatenate((w_all['w1m_bin0'], w_all['w1m_bin1'], w_all['w1m_bin2'], w_all['w1m_bin3']))
w2p_all = np.concatenate((w_all['w2p_bin0'], w_all['w2p_bin1'], w_all['w2p_bin2'], w_all['w2p_bin3']))
w2m_all = np.concatenate((w_all['w2m_bin0'], w_all['w2m_bin1'], w_all['w2m_bin2'], w_all['w2m_bin3']))

g1p_all_mean = np.average(g1p_all, weights=w1p_all)
g1m_all_mean = np.average(g1m_all, weights=w1m_all)
g2p_all_mean = np.average(g2p_all, weights=w2p_all)
g2m_all_mean = np.average(g2m_all, weights=w2m_all)

R12_all = g1p_all_mean - g1m_all_mean / 0.02
R21_all = g2p_all_mean - g2m_all_mean / 0.02

response_components_all = np.array([R12_all, R21_all])
np.savetxt(path+'/Response_non-diag_allbins_3oct.txt', response_components_all)