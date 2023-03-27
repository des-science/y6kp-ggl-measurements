import yaml
from destest import destest
import treecorr


# basic dict props
destest_dict_ = {
                 'output_exists' : True,
                 'use_mpi'       : False,
                 'source'        : 'hdf5',
                 'dg'            : 0.01
                }

# populates a full destest yaml dict for each catalog selection based on the limited catalog input info provided in the common cats.yaml file
def create_destest_yaml(params, name, cal_type, group, table, select_path):
    """
    Creates the input dictionary structure from a passed dictionary rather than reading from a yaml file
    """
    destest_dict                = destest_dict_.copy()
    destest_dict['load_cache']  = params['load_cache']
    destest_dict['output']      = params['output']
    destest_dict['name']        = name
    destest_dict['filename']    = params['datafile']
    destest_dict['param_file']  = params['param_file']
    destest_dict['cal_type']    = cal_type
    destest_dict['group']       = group
    destest_dict['table']       = table
    destest_dict['select_path'] = select_path
    destest_dict['e']           = ['e_1','e_2']
    destest_dict['Rg']          = ['R11','R22']
    destest_dict['w']           = 'weight'

    return destest_dict

# build selector (and calibrator) classes from destest for the catalog.
def load_catalog(pipe_params, name, cal_type, group, table, select_path, inherit=None, return_calibrator=None):
    """
    Loads data access and calibration classes from destest for a given yaml setup file
    """
    # input yaml file defining catalog
    params = create_destest_yaml(pipe_params, name, cal_type, group, table, select_path)

    # load destest source class to manage access to file
    source = destest.H5Source(params)

    # load destest selector class to manage access to data in a structured way
    if inherit is None:
        sel = destest.Selector(params,source)
    else:
        sel = destest.Selector(params,source,inherit=inherit)

    # load destest calibrator class to manage calibration of the catalog
    if return_calibrator is not None:
        cal = return_calibrator(params,sel)
        return sel, cal
    else:
        return sel