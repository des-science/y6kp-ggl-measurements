import sys
import yaml
sys.path.append('../../destest/')
import destest
import numpy as np

param_file = './destest_metacal.yaml'
params_mcal = yaml.load(open(param_file))
source_mcal = destest.H5Source(params_mcal)
selector_mcal = destest.Selector(params_mcal,source_mcal)

param_file = './destest_bpz.yaml'
params_bpz = yaml.load(open(param_file))
source_bpz = destest.H5Source(params_bpz)
selector_bpz = destest.Selector(params_bpz,source_bpz,inherit=selector_mcal)

param_file = './destest_gold.yaml'
params_gold = yaml.load(open(param_file))
source_gold = destest.H5Source(params_gold)
selector_gold  = destest.Selector(params_gold,source_gold,inherit=selector_mcal)

ra = selector_gold.get_col('ra') 
dec = selector_gold.get_col('dec') 
e1 = selector_mcal.get_col('e1')
e2 = selector_mcal.get_col('e2')
psf_e1 = selector_mcal.get_col('psf_e1')
psf_e2 = selector_mcal.get_col('psf_e2')
snr = selector_mcal.get_col('snr')
size = selector_mcal.get_col('size')

bpz_mean = selector_bpz.get_col('zmean_sof')
bpz_zmc = selector_bpz.get_col('zmc_sof')

calibrator = destest.MetaCalib(params_mcal,selector_mcal)
mask = selector_mcal.get_mask()
R11, _, _ = calibrator.calibrate('e1')
R22, _, _ = calibrator.calibrate('e2')
R = np.mean([R11, R22])
photoz_mask = (bpz_mean<0.43)&(bpz_mean>0.2)
R11, _, _ = calibrator.calibrate('e1', mask = photoz_mask)
