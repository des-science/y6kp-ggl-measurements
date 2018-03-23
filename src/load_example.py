import destest
import yaml
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

ra = selector_gold.get_col('alphawin_j2000') 
dec = selector_gold.get_col('deltawin_j2000') 
e1 = selector_mcal.get_col('e1')
e2 = selector_mcal.get_col('e2')
psf_e1 = selector_mcal.get_col('psf_e1')
psf_e2 = selector_mcal.get_col('psf_e2')
snr = selector_mcal.get_col('snr')
psf_size = selector_mcal.get_col('psf_size')
size_ratio = selector_mcal.get_col('size_ratio')

bpz_mean = selector_bpz.get_col('bpz_zmean_sof')
bpz_zmc = selector_bpz.get_col('bpz_zmc_sof')



