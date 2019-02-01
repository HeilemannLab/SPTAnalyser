# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 15:32:28 2019

@author: pcoffice37

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import time
import numpy as np
from pySPT.analysis import cell

def main():

    one_cell = cell.Cell()
    file_type = "palmtracer"
    if file_type == "swift":
        one_cell.pixel_size = 1
    elif file_type == "palmtracer":
        one_cell.pixel_size = 158
    start = time.time()
    one_cell.load_file()
    #cell.get_trajectories()
    one_cell.analyse_trajectories()
    end = time.time()
    print("Evalutation took {} seconds".format(end-start))
    #cell.plot_trajectorie(1)
    
    
if __name__ == "__main__":
    main()
    