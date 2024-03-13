"""
------------------------
Run the GGL measurements
------------------------
"""

import sys
import os

if __name__ == "__main__":
    
    # folder containing the input file with parameters
    if len(sys.argv) > 1:
        input_dir = sys.argv[1]
    else:
        input_dir = input_dir = './input'
    sys.path.append(input_dir)

    # code to run the measurements
    import ggly6 as ggl_run
    run = ggl_run.GGL(input_dir=input_dir)
    run.run_gammat_measurement()
    
    exit()
