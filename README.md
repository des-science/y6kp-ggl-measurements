# Galaxy-galaxy lensing measurements and tests
Credit: Giulia Giannini, Georgios Zacharegkas, Elisa Legnani ... (add your name here)

DES Y6 pipeline to validate galaxy-galaxy lensing datavectors.

### Install dependencies
See `requirements.txt`. They can all be installed at once by running:
```
pip install -r requirements.txt
```

### Run 
To execute, run:
```
pithon run_ggl.py input_params.yaml
```
Make sure you specify settings and file paths in your configuration file `input_params.yaml`. See the folder `input` for examples.

See the notebook `examples_run.ipynb` for examples of simple runs.

When running on the full data and including all the correction factors, the code will produce a 2pt file which requires blinding. 
The blinding repo is in https://github.com/des-science/y6-3x2pt/tree/main/blinding

### Read outputs
See the notebook `read_ggl_outputs.ipynb`.

------
Blinding directory:
Directory including ini files to use for the shift to apply to DV.

Parameters shifted:
- sigma8
- w

Source file:
d_l-sig8-hm20-tatt_y3fidCuts.test.sh

Example on how to run (specifically on nersc):
The following example scripts have the fiducial string keys for each of the two WL samples
- blinding_Y6_example_BFD.sh
- blinding_Y6_example_MDET.sh

Initialise the CSL_DIR variable to point to your cosmosis-standard-library

Make sure you include ANGLEMIN ANGLEMAX in your 2pt file

------
Please contact Giulia Giannini (giulia.giannini15@gmail.com) or Georgios Zacharegkas (gzacharegkas@uchicago.edu) of Elisa Legnani (elegnani@ifae.es) for comments or questions. 

