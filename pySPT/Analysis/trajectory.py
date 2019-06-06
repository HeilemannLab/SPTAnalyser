# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 16:14:28 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Create a trajectory object with informations of diffusion, diffusion type, length of trajectory, visualizations ...
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.stats import chisquare
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec

class Trajectory():
    def __init__(self, locs, tau_thresh, camera_dt, degree, min_D, points_D):
        self.trajectory_number = 0  # equals first column of localization
        self.MSDs = []  # stores all MSD values
        self.times = []  # stores all time steps for MSD values
        self.MSD_fit = []  # stores values of fit area of times, MSDs, fit and residues
        self.MSD_D = []  # col0 = times, col1 = MSD vals, col2 = fit, col3= res for the first 4 values
        self.localizations = locs  # np.array with col0 = track id, col1 = frames, col2 = x [nm], col3 = y [nm], col4 = placeholder, col5 = intensity, col6 = seg id
        self.dt = camera_dt  # camera integration time in s
        self.dof = degree  # degree of freedom (to determin D)
        self.D_min = min_D  # minimum of diffusion coeff to be detected [micrometer^2/s]
        self.length_MSD = 0
        self.length_trajectory = 0
        self.D = 0.0  # diffusion coefficient
        self.dD = 0.0  # error of diffusion coefficient
        self.chi_D = 0.0  # chi^2 of linear fit for MSD
        self.chi_MSD_fit = 0.0  # chi^2 for the MSD 60% fit
        self.MSD_0 = 0.0  # y-intercept
        self.dMSD_0 = 0.0  # not determined yet ...
        self.fit_area = 0.6  # 60 % of the MSD plot will be fitted
        self.tau = 0.0  # tau value, derived by rossier fit parameters Dconfined and r
        self.dtau = 10.0  # error of tau value, derived with gauß
        self.D_conf = 0.01  # confined diffusion
        self.dD_conf = 0.0  # error of confined diffusion, determined by covariance matrix of fit
        self.r = 0.01  # confinement radius
        self.dr = 0.0  # error of confinement radius, determined by covariance matrix of fit
        self.tau_threshold = tau_thresh  # check if free or confined: tau > tau_thr -> free, tau < tau_thr -> confined
        self.immobility = False
        self.confined = True  # if confinded false -> free true
        self.analyse_successful = True
        self.points_fit_D = points_D  # number of points that will be fittet to extract D
        
    def calc_trajectory_number(self):
        self.trajectory_number = self.localizations[0,0]   
        
    def calc_length_trajectory(self):
        self.length_trajectory = np.shape(self.localizations)[0]
        
    def calc_length_MSD(self):
        self.length_MSD = np.shape(self.localizations)[0] - 1
        
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
                MSD = (self.localizations[i+time,2] - self.localizations[i,2])**2 + (self.localizations[i+time,3] - self.localizations[i,3])**2
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
        self.MSD_D = np.zeros([self.points_fit_D,4])
        self.MSD_D[:,1] = self.MSDs[:self.points_fit_D]
        for i in range(self.points_fit_D):
            self.MSD_D[i,0] = (i+1)*self.dt
        slope, intercept, r_value, p_value, std_err = linregress(self.MSD_D[:,0], self.MSD_D[:,1])
        self.MSD_D[:,2] = self.function_linear_fit(self.MSD_D[:,0], slope, intercept)
        self.MSD_D[:,3] = self.MSD_D[:,1] - self.MSD_D[:,2]
        self.D = slope/self.dof
        [chisq, p] = chisquare(self.MSD_D[:,1], self.function_linear_fit(self.MSD_D[:,0], slope, intercept))
        self.chi_D = chisq
# =============================================================================
#         if not (self.D > 1.0*10**(-5.0)):
#             self.D = 1.0*10**(-5.0)
# =============================================================================
        self.dD = std_err/self.dof
        self.MSD_0 = intercept
        self.dMSD_0 = ""
        
    def function_linear_fit(self, times, slope, intercept):
        return times*slope + intercept    
        
    def plot_diffusion(self):
        #x1, x2 = self.MSD_D[:,0].min(), self.MSD_D[:,0].max()  # x1 = min, x2 = max
        #sp1_y1, sp1_y2 = self.MSD_D[:,1].min(), self.MSD_D[:,1].max()
        #sp2_y1, sp2_y2 = self.MSD_D[:,3].min(), self.MSD_D[:,3].max()
        #fig = plt.figure()
        gridspec.GridSpec(4,4)  # set up supbplot grid 
        sp_1 = plt.subplot2grid((4,4), (0,0), colspan=4, rowspan=3)  # start left top = (0,0) = (row,column)
        sp_1.tick_params(axis='x',  # changes apply to the x-axis
                        which='both',  # both major and minor ticks are affected
                        bottom=False,  # ticks along the bottom edge are off
                        top=False,  # ticks along the top edge are off
                        labelbottom=False) # labels along the bottom edge are off
        #sp_1 = fig.add_subplot(2, 1, 1)  # (row, column, index)
        sp_1.plot(self.MSD_D[:,0], self.MSD_D[:,1], "o", color="0.5", label="MSD values")
        sp_1.plot(self.MSD_D[:,0], self.MSD_D[:,2], "--c", label="linear fit")
        sp_1.legend()
        sp_1.set_title("MSD-Plot: Diffusion")
        #sp_1.set_xlabel("Number of data points used in MJD calculation")
        sp_1.set_ylabel("MSD")
        #sp_1.axis((x1, x2, sp1_y1, sp1_y2))
        sp_2 = plt.subplot2grid((4,4), (3,0), colspan=4, rowspan=1)
        #sp_2 = fig.add_subplot(2, 1, 2) 
        residue_line = np.zeros(len(self.MSD_D[:,0]))
        sp_2.plot(self.MSD_D[:,0], residue_line, ":", color = "0.75")
        sp_2.plot(self.MSD_D[:,0], self.MSD_D[:,3], "*", color="0.5", label="residues")
        sp_2.legend()
        sp_2.set_ylabel("Residue")
        sp_2.set_xlabel("Time step [s]")  # Number of data points used in MJD calculation
        #sp_2.axis((x1, x2, sp2_y1, sp2_y2))
        plt.show() 
        
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
            
    def plot_full_MSD_immob(self):
        """
        Plot 60% of MSD if molecule is immob without rossier fit.
        """
        self.create_MSD_values()
        #x1, x2 = 0, self.MSD_fit[:,0].max()  # x1 = min, x2 = max
        #y1, y2 = 0, self.MSD_fit[:,1].max()
        plt.plot(self.MSD_fit[:,0], self.MSD_fit[:,1], "o", color = "0.5", label="MSD values")
        plt.legend()
        plt.title("MSD-Plot (60 % of values)")
        plt.xlabel("Time step [s]")
        plt.ylabel("MSD")
        plt.show()
        #plt.axis((x1,x2,y1,y2))
        
    def create_MSD_values(self):
        max_times = np.rint(self.length_MSD*self.fit_area)      
        times = np.arange(1, max_times+1, 1.0)
        times[:] = times[:] * self.dt
        MSDs = self.MSDs[:int(max_times)]
        self.MSD_fit = np.zeros([len(times),4])
        self.MSD_60 = np.zeros([len(times),2])
        self.MSD_fit[:,0] = times
        self.MSD_fit[:,1] = MSDs     
        
# =============================================================================
#     def function_full_MSD_tau(self, t, r, tau):
#         return (4.0*r**2.0)/3.0*(1.0-np.exp(-t/tau))
#         
#     def fit_full_MSD_tau(self):
#         max_times = np.rint(self.length_MSD*self.fit_area)      
#         times = np.arange(1, max_times+1, 1.0)
#         times[:] = times[:] * self.dt
#         MSDs = self.MSDs[:int(max_times)]
#         self.MSD_fit = np.zeros([len(times),4])
#         self.MSD_60 = np.zeros([len(times),2])
#         self.MSD_fit[:,0] = times
#         self.MSD_fit[:,1] = MSDs
#         try:
#             [self.r, self.tau], self.cov = curve_fit(self.function_full_MSD, times, MSDs, p0 = (self.r, self.tau), method = "lm")
#         except:
#             self.analyse_successful = False
#         if self.analyse_successful:
#             [chisq, p] = chisquare(MSDs, self.function_full_MSD(times, self.r, self.tau))
#             self.chi_MSD_fit = chisq
#             self.D_conf = self.r**2.0/(3.0*self.tau)
#             self.dr = self.cov[0,0]
#             self.dtau = self.cov[1,1]
#             self.MSD_fit[:,2] = self.function_full_MSD(times, self.r, self.tau)
#             self.MSD_fit[:,3] = self.MSD_fit[:,1] - self.MSD_fit[:,2]
# =============================================================================
            
    def fit_full_MSD(self):
        max_times = np.rint(self.length_MSD*self.fit_area)      
        times = np.arange(1, max_times+1, 1.0)
        times[:] = times[:] * self.dt
        MSDs = self.MSDs[:int(max_times)]
        self.MSD_fit = np.zeros([len(times),4])
        self.MSD_60 = np.zeros([len(times),2])
        self.MSD_fit[:,0] = times
        self.MSD_fit[:,1] = MSDs
        try:
            [self.r, self.D_conf], self.cov = curve_fit(self.function_full_MSD, times, MSDs, p0 = (self.r, self.D_conf), method = "lm")
        except:
            self.analyse_successful = False
        if self.analyse_successful:
            [chisq, p] = chisquare(MSDs, self.function_full_MSD(times, self.r, self.D_conf))
            self.chi_MSD_fit = chisq
            self.tau = self.r**2/(3*self.D_conf)
            self.dr = self.cov[0,0]
            self.dD_conf = self.cov[1,1]  
            self.dtau = math.sqrt((2*self.r/(3*self.D_conf)*self.dr)**2 + (-self.r**2/(3*self.D_conf**2)*self.dD_conf)**2)
            self.MSD_fit[:,2] = self.function_full_MSD(times, self.r, self.D_conf)
            self.MSD_fit[:,3] = self.MSD_fit[:,1] - self.MSD_fit[:,2]
            if self.dr == math.inf:
                self.analyse_successful = False
            
    def function_full_MSD(self, t, r, D):
        return (4.0*r**2.0)/3.0*(1.0-np.exp(-t*3*D/r**2.0))
        
    def plot_full_MSD(self):
        #x1, x2 = 0, self.MSD_fit[:,0].max()  # x1 = min, x2 = max
        #sp1_y1, sp1_y2 = 0, self.MSD_fit[:,1].max()
        #sp2_y1, sp2_y2 = self.MSD_fit[:,3].min(), self.MSD_fit[:,3].max()
        #fig = plt.figure()
        gridspec.GridSpec(4,4)  # set up supbplot grid 
        sp_1 = plt.subplot2grid((4,4), (0,0), colspan=4, rowspan=3)  # start left top = (0,0) = (row,column)
        sp_1.tick_params(axis='x',  # changes apply to the x-axis
                        which='both',  # both major and minor ticks are affected
                        bottom=False,  # ticks along the bottom edge are off
                        top=False,  # ticks along the top edge are off
                        labelbottom=False) # labels along the bottom edge are off
        #sp_1 = fig.add_subplot(2, 1, 1)  # (row, column, index)
        sp_1.plot(self.MSD_fit[:,0], self.MSD_fit[:,1], "o", color="0.5", label="MSD values")
        sp_1.plot(self.MSD_fit[:,0], self.MSD_fit[:,2], "--c", label="rossier fit")
        #sp_1.plot(self.MSD_fit_ML[:,0], self.MSD_fit_ML[:,2], "--m", label="origin fit")
        sp_1.legend()
        sp_1.set_title("MSD-Plot: Type")
        #sp_1.set_xlabel("Number of data points used in MJD calculation")
        sp_1.set_ylabel("MSD")
        #sp_1.axis((x1, x2, sp1_y1, sp1_y2))
        sp_2 = plt.subplot2grid((4,4), (3,0), colspan=4, rowspan=1)
        #sp_2 = fig.add_subplot(2, 1, 2) 
        residue_line = np.zeros(len(self.MSD_fit[:,0]))
        sp_2.plot(self.MSD_fit[:,0], residue_line, ":", color = "0.75")
        sp_2.plot(self.MSD_fit[:,0], self.MSD_fit[:,3], "*", color = "0.5", label="residues")
        sp_2.legend()
        sp_2.set_ylabel("Residue")
        sp_2.set_xlabel("Time step [s]")  # Number of data points used in MJD calculation
        #sp_2.axis((x1, x2, sp2_y1, sp2_y2))
        plt.show() 
        
    def show_trajectory(self):
        plt.plot(self.localizations[:,2], self.localizations[:,3], linestyle="--", marker="o", label="localization")
        plt.legend()
        plt.title("Trajectory of one particle")
        plt.xlabel(r"$\mu$" + "m in x direction")
        plt.ylabel(r"$\mu$" + "m in y direction")
        plt.show()

    def analyse_particle(self):
        self.calc_trajectory_number()
        self.calc_length_trajectory()
        self.calc_length_MSD()
        self.calc_MSD()
        self.calc_diffusion()
        self.check_immobility()
        self.create_MSD_values()
        if not (self.immobility):
            self.fit_full_MSD()
            self.check_confined()

    def print_particle(self):
        print("Number:", int(self.trajectory_number))
        print("Trajectory length:", int(self.length_trajectory))
        print("Diffusion coefficient: {} \u03BCm\u00b2/s".format(self.D))
        print("MSD0: {} \u03BCm\u00b2".format(self.MSD_0))
        print("chi\u00b2 linear fit: {} \u03BCm\u2074".format(self.chi_D))
        print("Type immobile:", self.immobility)
        if not self.immobility:
            print("Analyse successful?", self.analyse_successful)
            if self.analyse_successful:
                print("chi\u00b2 rossier fit: {} \u03BCm\u2074".format(self.chi_MSD_fit))
                print("Type confined:", self.confined)
                print("Type free:", not self.confined)
                print("D_conf: {} \u03BCm\u00b2/s".format(self.D_conf))
                print("r_conf: {} \u03BCm".format(self.r))
                print("tau: {} s".format(self.tau))
                print("tau threshold: {} s".format(self.tau_threshold))

    def plot_particle(self):
        self.show_trajectory()
        self.plot_diffusion()
        if not self.immobility:
            self.plot_full_MSD()
        elif self.immobility:
            self.plot_full_MSD_immob()
        self.print_particle()
   
   
def main():
    trajectory = Trajectory()
    trajectory.analyse_particle()
    trajectory.plot_particle()
        
    
if __name__ == "__main__":
    main()
    