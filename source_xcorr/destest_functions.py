import yaml
import sys
#sys.path.append('../../destest/')
#import destest
from destest import destest
import treecorr
from info import config


# basic dict props
destest_dict_ = {
    'output_exists' : True,
    'use_mpi'       : False,
    'source'        : 'hdf5',
    'dg'            : 0.01
    }

# Populates a full destest yaml dict for each catalog selection based on the limited catalog input info provided in the common cats.yaml file
def create_destest_yaml( params, name, cal_type, group, table, select_path ):
    """
    Creates the input dictionary structure from a passed dictionary rather than reading froma yaml file.
    """

    destest_dict = destest_dict_.copy()
    destest_dict['load_cache'] = params['load_cache']
    destest_dict['output'] = params['output']
    destest_dict['name'] = name
    destest_dict['filename'] = config['filename_mastercat'] # from info.py to make sure we are consistent
    destest_dict['param_file'] = params['param_file']
    destest_dict['cal_type'] = cal_type
    destest_dict['group'] = group
    destest_dict['table'] = table
    destest_dict['select_path'] = select_path
    destest_dict['e'] = ['e_1','e_2']
    destest_dict['Rg'] = ['R11','R22']
    destest_dict['w'] = 'weight'

    return destest_dict

# Build selector (and calibrator) classes from destest for the catalog.
def load_catalog(pipe_params, name, cal_type, group, table, select_path, inherit=None, return_calibrator=None):
    """
    Loads data access and calibration classes from destest for a given yaml setup file.
    """

    # Input yaml file defining catalog
    params = create_destest_yaml(pipe_params, name, cal_type, group, table, select_path)

    # Load destest source class to manage access to file
    source = destest.H5Source(params)

    # Load destest selector class to manage access to data in a structured way
    if inherit is None:
        sel = destest.Selector(params,source)
    else:
        sel = destest.Selector(params,source,inherit=inherit)

    # Load destest calibrator class to manage calibration of the catalog
    if return_calibrator is not None:
        cal = return_calibrator(params,sel)
        return sel, cal
    else:
        return sel

'''
# Read yaml file that defines all the catalog selections used
params = yaml.load(open('cats.yaml'))
params['param_file'] = 'cats.yaml'

# Source catalog
source_selector, source_calibrator = load_catalog(
    params, 'mcal', 'mcal', params['source_group'], params['source_table'], params['source_path'], return_calibrator=destest.MetaCalib)

# Redmagic catalog
lens_selector, lens_calibrator = load_catalog(
    params, 'rm', None, params['lens_group'], params['lens_table'], params['lens_path'], return_calibrator=destest.NoCalib)

# Mag lim catalog
alt_lens_selector, alt_lens_calibrator = load_catalog(
    params, 'maglim', None, params['alt_lens_group'], params['alt_lens_table'], params['alt_lens_path'], return_calibrator=destest.NoCalib)

# Gold catalog
gold_selector = load_catalog(
    params, 'gold', 'mcal', params['gold_group'], params['gold_table'], params['gold_path'], inherit=source_selector)

# BPZ (or DNF) catalog, depending on paths in cats.yaml file (exchange bpz and dnf)
pz_selector = load_catalog(
    params, 'pz', 'mcal', params['pz_group'], params['pz_table'], params['pz_path'], inherit=source_selector)

# Redmagic random catalog
ran_selector = load_catalog(
    params, 'ran', None, params['ran_group'], params['ran_table'], params['ran_path'])

# Get some source photo-z binning information, cut to range 0.1<z_mean<1.3                                                  
for i in range(4):
    pzbin = pz_selector.get_col('bhat') # 5-tuple for metacal (un)sheared versions                                         
    mask = [pzbin[j] == i for j in range(5)] # First tomographic bin                                                               
    # Note that get_col() returns a tuple. If its a catalog like gold, it will have length 0, but for something like metacal, it will have length 5 (in the order of the table variable list passed in cats.yaml, i.e., 'unsheared', 'sheared_1p', 'sheared_1m', 'sheared_2p', 'sheared_2m')                                                               
    # Note that get_col() applies the index mask specified by the 'path' variable in the cats.yaml file automatically.         # Get responses (c, which doesn't exist for our catalogs), and weights                                                      
    R1,c,w = source_calibrator.calibrate('e_1', mask=mask) # Optionally pass an additional mask to use when calculating the selection response. The returned R1 is <Rg_1 + Rs_1>. To get an array of R's, use return_wRg=True to get [Rg_1+Rg_2]/2 for each object or return_wRgS=True to include the selection response. return_full=True returns the non-component-averaged version of the full response.
    print(R1,c,w)
    g1=source_selector.get_col('e_1')
    print(len(w),[len(g1[j]) for j in range(4)])
    R2,c,w = source_calibrator.calibrate('e_2', mask=mask)
    print(R2,c,w)

# Load ra,dec from gold catalog
ra  = gold_selector.get_col('ra')[0]
dec = gold_selector.get_col('dec')[0]

# Load ra,dec from redmagic catalog
ra  = lens_selector.get_col('ra')[0]
dec = lens_selector.get_col('dec')[0]
ra2  = ran_selector.get_col('ra')[0]
dec2 = ran_selector.get_col('dec')[0]

# Get e1,e2 (all five version)
g1=source_selector.get_col('e_1')
g2=source_selector.get_col('e_2')

pz_selector = load_catalog(
        params, 'pzdnf', None, params['dnf_group'], params['dnf_table'], params['dnf_path'], inherit=alt_lens_selector)

'''
