# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 11:03:55 2019

@author: pcoffice37

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.stats import chisquare
from scipy.optimize import curve_fit

class Molecule():
    def __init__(self):
        self.MSDs = []  # stores all MSD values
        self.times = []  # stores time steps for MSD values
        self.localizations = []  # list with tuples of localizations (x1,y1), (x2,y2) ...
        self.array_trajectory = []
        self.array_MSD = []
        self.sorted_molecules = [] # col0 = molecule, col1 = frames, col2 = x, col3 = y, col4 = intensity
        self.dt = 0.02  # camera integration time in s
        self.pixel_size = 158  # [nm] multiply with palmtracer, swift = 1 -> in nm -> *10^-3 micrometer!
        self.dof = 4  # degree of freedom 
        self.D_min = 0.0065  # [micrometer^2/s]
        self.length_MSD = 0
        self.length_trajectory = 0
        self.D = 0.0  # diffusion coefficient
        self.dD = 0.0
        self.chi = 0.0  # chi^2 of linear fit for MSD
        self.MSD_0 = 0.0  # y-intercept
        self.dMSD_0 = 0.0
        self.tau = 0.01  # fit for free confined
        self.dtau = 10.0
        self.D_conf = 0.01  # confinded diffusion
        self.r = 0.01  # confinement radius
        self.dr = 0.0
        self.tau_threshold = 0.12  # check if free or confined
        self.immobility = False
        self.confined = True  # if confinded false -> free true
        self.analyse_succesfull = True
        
    def load_file(self):
        #file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\sorted.txt"
        file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\cell_1_MMStack_Pos0.ome.MIA\\tracking\\cell_1_MMStack_Pos0.ome_MIA.trc"
        self.sorted_molecules = np.loadtxt(file_name, usecols = (0, 1, 2, 3, 5)) # col0 = molecule, col1 = frames, col2 = x, col3 = y, col4 = intensity
        self.sorted_molecules[:,2] = self.sorted_molecules[:,2] * self.pixel_size *10**-3
        self.sorted_molecules[:,3] = self.sorted_molecules[:,3] * self.pixel_size *10**-3
        #print(self.sorted_molecules)
        
    def get_localization(self):
        # get localizations for first molecule as xy tuples in list
        self.localizations = []  # (x1, y1), (x2, y2) ...
        #print(self.sorted_molecules[20,0])
        for i in range(0, len(self.sorted_molecules[:,0])):
            if self.sorted_molecules[i,0] == 206:
                self.localizations.append((self.sorted_molecules[i,2], self.sorted_molecules[i,3]))
        
    def calc_length_trajectory(self):
        self.length_trajectory = len(self.localizations)
        
    def calc_length_MSD(self):
        self.length_MSD = len(self.localizations) - 1
        
    def calc_MSD(self):
        """
        Creat array with mean MSD-values for different time intervals dt.
        Calculate all possible MSD-values for certain time -> MSD_dt.
        The mean will be stored in MSDs array.
        """
        self.MSDs = np.zeros(self.length_MSD)
        for time in range(1, self.length_trajectory):
            MSD_time = np.zeros(self.length_trajectory - time)
            for i in range(0, self.length_trajectory - time):
                MSD = (self.localizations[i+time][0] - self.localizations[i][0])**2 + (self.localizations[i+time][1] - self.localizations[i][1])**2
                MSD_time[i] = MSD
            self.MSDs[time-1] = MSD_time.mean()
        self.times = np.arange(1, self.length_MSD+1, 1.0)
        self.times[:] = self.times[:] * self.dt
            
    def calc_diffusion(self):
        """
        Calculate diffusion coefficient based on first 4 MSD values
        times = first 4 time distances for MSD calculation -> *dt
        if diffusion is 2D -> D = slope/4, if 3D D = slope/6 ...
        """
        MSDs = self.MSDs[:4]
        times = [1*self.dt, 2*self.dt, 3*self.dt, 4*self.dt]
        slope, intercept, r_value, p_value, std_err = linregress(times, MSDs)
        self.D = slope/self.dof
        if not (self.D > 1.0*10**(-5.0)):  
            self.D = 1.0*10**(-5.0)
        self.dD = std_err/self.dof
        self.MSD_0 = intercept
        self.dMSD_0 = ""
        fit_fn = np.poly1d([slope, intercept])
        plt.title("MSD-Plot")
        plt.ylabel("MSD")
        plt.xlabel("Time step [s]")
        plt.plot(times, MSDs, "yo", label="MSD values")
        plt.plot(times, fit_fn(times), "--k", label="linear fit")
        plt.legend()
        [chisq, p] = chisquare(MSDs, fit_fn(times))
        self.chi = chisq
        plt.show()
        print("Diffusion coeff:", self.D)
        print("MSD0:", self.MSD_0)
 
    def check_immobility(self):
        self.immobility = False
        if self.D < self.D_min:
            self.immobility = True     

    def check_confined(self):
        """
        Half time interval is the rounded up time at 30 % of the full time spawn.
        Confined if tau < tau threshold -> return True for self.confined
        Free if tau > tau threshold -> return False
        """
        self.confined = False
        if self.tau < self.tau_threshold:
            self.confined = True
    
    def function_full_MSD(self, t, r, tau):
        return (4.0*r**2.0)/3.0*(1.0-np.exp(-t/tau))
    
    def fit_full_MSD(self):
        max_times = np.ceil(self.length_MSD*0.6)      
        times = np.arange(1, max_times+1, 1.0)
        times[:] = times[:] * self.dt
        MSDs = self.MSDs[:int(max_times)]
        #MSDs = self.MSDs[:max_times]
        print("Full length MSD:", self.length_MSD)
        print("Fitted length MSD (60 %):", len(MSDs))
        try:
            [self.r, self.tau], self.cov = curve_fit(self.function_full_MSD, times, MSDs, p0 = (self.r, self.tau), method = "lm")
        except:
            self.analyse_succesfull = False
        if self.analyse_succesfull:
            self.D_conf = self.r**2.0/(3.0*self.tau)
            self.dr = self.cov[0,0]
            self.dtau = self.cov[1,1]
            #[chisq, p] = chisquare(MSDs, fit_fn(times))
        print("Analysis succes?", self.analyse_succesfull)        
        
    def plot_full_MSD(self):
        max_times = np.ceil(self.length_MSD*0.6)
        times = np.arange(1, max_times+1, 1.0)
        times[:] = times[:] * self.dt
        MSDs = self.MSDs[:int(max_times)]
        MSD_fit = self.function_full_MSD(times, self.r, self.tau)
        MSD_fit_ML = self.function_full_MSD(times, 0.05716, 0.000385009101495349)
        fig = plt.figure()
        sp = fig.add_subplot(1,1,1)
        sp.plot(times, MSDs, "o")
        sp.plot(times, MSD_fit, "--k")
        sp.plot(times, MSD_fit_ML, "--g")
        plt.show()
        
    def analyse_particle(self):
        self.calc_length_trajectory()
        self.calc_length_MSD()
        self.calc_MSD()
        self.calc_diffusion()
        self.check_immobility()
        if not (self.immobility):
            self.fit_full_MSD()
            self.check_confined()
        self.plot_full_MSD()
        print("Type immobile:", self.immobility)
        print("Type confined:", self.confined)
        print("rconf:", self.r)
        print("tau:", self.tau)
   

def main():
    molecule = Molecule()
    molecule.load_file()
    molecule.get_localization()
    molecule.analyse_particle()
    
    
if __name__ == "__main__":
    main()
    