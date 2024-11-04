import numpy as np
import pandas as pd

#Load all ellipticity files

g = {}
path = "/global/cfs/projectdirs/des/awhyley/ggl/response_tests/"

for zbin in range(4):
    print("loading z bin")
    g['g1p_bin{}'.format(zbin)] = np.loadtxt(path+'/g1p_bin{}.txt'.format(zbin))
    g['g1m_bin{}'.format(zbin)] = np.loadtxt(path+'/g1m_bin{}.txt'.format(zbin))
    g['g2p_bin{}'.format(zbin)] = np.loadtxt(path+'/g2p_bin{}.txt'.format(zbin))
    g['g2m_bin{}'.format(zbin)] = np.loadtxt(path+'/g2m_bin{}.txt'.format(zbin))
    
#Average over all z bins

print("averaging ellipticities over all bins")

g1p_all = np.concatenate((g['g1p_bin0'], g['g1p_bin1'], g['g1p_bin2'], g['g1p_bin3']))
g1m_all = np.concatenate((g['g1m_bin0'], g['g1m_bin1'], g['g1m_bin2'], g['g1m_bin3']))
g2p_all = np.concatenate((g['g2p_bin0'], g['g2p_bin1'], g['g2p_bin2'], g['g2p_bin3']))
g2m_all = np.concatenate((g['g2m_bin0'], g['g2m_bin1'], g['g2m_bin2'], g['g2m_bin3']))

g1p_all_mean = np.mean(g1p_all)
g1m_all_mean = np.mean(g1m_all)
g2p_all_mean = np.mean(g2p_all)
g2m_all_mean = np.mean(g2m_all)

print("means:")
print(g1p_all_mean, g1m_all_mean, g2p_all_mean, g2m_all_mean)

#Calculate response matrix components

