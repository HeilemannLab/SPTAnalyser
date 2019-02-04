# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 11:45:55 2019

@author: pcoffice37

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

from . import trajectory
from . import cell

class TrajectoryStatistics():
    def __init__(self):
        self.trajectory_list = []
        
    def plot_trajectory(self, cell, number):
        cell = int(cell) - 1
        number = int(number) -1
        self.trajectory_list[cell][number].plot_particle()
        
    
        
    
    
    
    
def main():
    pass
        
    
if __name__ == "__main__":
    main()
        
    
    