import sys
import yaml
sys.path.append('../../destest/')
import destest

param_file = './destest_redmagic.yaml'
print param_file
params_rm = yaml.load(open(param_file))
source_rm = destest.H5Source(params_rm)
selector_rm = destest.Selector(params_rm,source_rm)


ra = selector_rm.get_col('ra') 
dec = selector_rm.get_col('dec') 
print len(ra)

