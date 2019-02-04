# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 16:14:28 2019

@author: pcoffice37

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""

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
import matplotlib.gridspec as gridspec
#from matplotlib import rc

class Trajectory():
    def __init__(self, locs):
        #self.trajectory_number = trajectory_number  # number of trajectory
        #self.len_all_trajectories = len_all_trajectories  # sum of length of all trajectories
        #self.all_trajectories = all_trajectories  # col0 = molecule, col1 = frames, col2 = x, col3 = y, col4 = intensity
        self.cell_number = 0  
        self.trajectory_number = 0  # equals first column of localization
        self.MSDs = []  # stores all MSD values
        self.times = []  # stores all time steps for MSD values
        self.MSD_fit = []  # stores values of fit area of times, MSDs, fit and residues

        self.MSD_D = []  # stores first 4 values for diffusion calc, times MSD, fit residues
        self.localizations = locs  # np.array with col0 = molecule, col1 = frames, col2 = x, col3 = y, col4 = intensity
        #self.sorted_molecules = [] # col0 = molecule, col1 = frames, col2 = x, col3 = y, col4 = intensity
        self.dt = 0.02  # camera integration time in s
        #self.pixel_size = 158  # [nm] multiply with palmtracer, swift = 1 -> in nm -> *10^-3 micrometer!
        self.dof = 4  # degree of freedom 
        self.D_min = 0.0065  # [micrometer^2/s]
        self.length_MSD = 0
        self.length_trajectory = 0
        self.D = 0.0  # diffusion coefficient
        self.dD = 0.0
        self.chi_D = 0.0  # chi^2 of linear fit for MSD
        self.chi_MSD_fit = 0.0  # chi^2 for the MSD 60% fit
        self.MSD_0 = 0.0  # y-intercept
        self.dMSD_0 = 0.0  # not determined yet ...
        self.fit_area = 0.6  # 60 % of the MSD plot will be fitted
        self.tau = 0.01  # fit for free confined
        self.dtau = 10.0
        self.D_conf = 0.01  # confinded diffusion
        self.dD_conf = 0.0
        self.r = 0.01  # confinement radius
        self.dr = 0.0
        self.tau_threshold = 0.12  # check if free or confined
        self.immobility = False
        self.confined = True  # if confinded false -> free true
        self.analyse_successful = True
        
        self.MSD_fit_ML= []  # comparing with origin fit
        self.D_conf_ML = 0.0
        self.tau_ML = 0.0
        self.r_ML = 0.0
        
# =============================================================================
#     def load_file(self):
#         #file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\sorted.txt"
#         file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\cell_1_MMStack_Pos0.ome.MIA\\tracking\\cell_1_MMStack_Pos0.ome_MIA.trc"
#         self.sorted_molecules = np.loadtxt(file_name, usecols = (0, 1, 2, 3, 5)) # col0 = molecule, col1 = frames, col2 = x, col3 = y, col4 = intensity
#         self.sorted_molecules[:,2] = self.sorted_molecules[:,2] * self.pixel_size *10**-3
#         self.sorted_molecules[:,3] = self.sorted_molecules[:,3] * self.pixel_size *10**-3
#         #print(self.sorted_molecules)
# =============================================================================
        
# =============================================================================
#     def get_localization(self):       
#         for i in range(0, self.len_all_trajectories):
#             if self.all_trajectories[i,0] == self.trajectory_number:
#                 self.localizations.append((self.all_trajectories[i,2], self.all_trajectories[i,3]))
# =============================================================================
# =============================================================================
#         print(self.localizations)
# =============================================================================
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
        self.MSD_D = np.zeros([4,4])
        self.MSD_D[:,1] = self.MSDs[:4]
        self.MSD_D[:,0] = [1*self.dt, 2*self.dt, 3*self.dt, 4*self.dt]
        slope, intercept, r_value, p_value, std_err = linregress(self.MSD_D[:,0], self.MSD_D[:,1])
        self.MSD_D[:,2] = self.function_linear_fit(self.MSD_D[:,0], slope, intercept)
        self.MSD_D[:,3] = self.MSD_D[:,1] - self.MSD_D[:,2]
        self.D = slope/self.dof
        [chisq, p] = chisquare(self.MSD_D[:,1], self.function_linear_fit(self.MSD_D[:,0], slope, intercept))
        self.chi_D = chisq
        if not (self.D > 1.0*10**(-5.0)):  
            self.D = 1.0*10**(-5.0)
        self.dD = std_err/self.dof
        self.MSD_0 = intercept
        self.dMSD_0 = ""
        
    def function_linear_fit(self, times, slope, intercept):
        return times*slope + intercept    
        
    def plot_diffusion(self):
        x1, x2 = self.MSD_D[:,0].min(), self.MSD_D[:,0].max()  # x1 = min, x2 = max
        sp1_y1, sp1_y2 = self.MSD_D[:,1].min(), self.MSD_D[:,1].max()
        sp2_y1, sp2_y2 = self.MSD_D[:,3].min(), self.MSD_D[:,3].max()
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
        sp_1.axis((x1, x2, sp1_y1, sp1_y2))
        sp_2 = plt.subplot2grid((4,4), (3,0), colspan=4, rowspan=1)
        #sp_2 = fig.add_subplot(2, 1, 2) 
        residue_line = np.zeros(len(self.MSD_D[:,0]))
        sp_2.plot(self.MSD_D[:,0], residue_line, ":", color = "0.75")
        sp_2.plot(self.MSD_D[:,0], self.MSD_D[:,3], "*", color="0.5", label="residues")
        sp_2.legend()
        sp_2.set_ylabel("Residue")
        sp_2.set_xlabel("Time step [s]")  # Number of data points used in MJD calculation
        sp_2.axis((x1, x2, sp2_y1, sp2_y2))
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
            
    def function_full_MSD_ML(self, t, r, D):
        return (4.0*r**2.0)/3.0*(1.0-np.exp(-t*3*D/r**2.0))
    
    def function_full_MSD(self, t, r, tau):
        return (4.0*r**2.0)/3.0*(1.0-np.exp(-t/tau))
    
    def fit_full_MSD_ML(self):
        max_times = np.rint(self.length_MSD*self.fit_area)      
        times = np.arange(1, max_times+1, 1.0)
        times[:] = times[:] * self.dt
        MSDs = self.MSDs[:int(max_times)]
        self.MSD_fit_ML = np.zeros([len(times),4])
        self.MSD_60 = np.zeros([len(times),2])
        self.MSD_fit_ML[:,0] = times
        self.MSD_fit_ML[:,1] = MSDs
        self.r_ML = 0.16891
        self.tau_ML = 0.01178
        self.D_conf_ML = 3*self.tau_ML/self.r_ML**2  # used wrong formular
        self.MSD_fit_ML[:,2] = self.function_full_MSD_ML(times, self.r_ML, self.D_conf_ML)
        self.MSD_fit_ML[:,3] = self.MSD_fit_ML[:,1] - self.MSD_fit_ML[:,2]
        
    def plot_full_MSD_immob(self):
        """
        Plot 60% of MSD if molecule is immob without rossier fit.
        """
        max_times = np.rint(self.length_MSD*self.fit_area)      
        times = np.arange(1, max_times+1, 1.0)
        times[:] = times[:] * self.dt
        MSDs = self.MSDs[:int(max_times)]
        self.MSD_fit = np.zeros([len(times),4])
        self.MSD_60 = np.zeros([len(times),2])
        self.MSD_fit[:,0] = times
        self.MSD_fit[:,1] = MSDs
        x1, x2 = 0, self.MSD_fit[:,0].max()  # x1 = min, x2 = max
        y1, y2 = 0, self.MSD_fit[:,1].max()
        plt.plot(self.MSD_fit[:,0], self.MSD_fit[:,1], "o", color = "0.5", label="MSD values")
        plt.legend()
        plt.title("MSD-Plot")
        plt.xlabel("Time step [s]")
        plt.ylabel("MSD")
        plt.axis((x1,x2,y1,y2))
        
    def fit_full_MSD(self):
        max_times = np.rint(self.length_MSD*self.fit_area)      
        times = np.arange(1, max_times+1, 1.0)
        times[:] = times[:] * self.dt
        MSDs = self.MSDs[:int(max_times)]
        self.MSD_fit = np.zeros([len(times),4])
        self.MSD_60 = np.zeros([len(times),2])
        self.MSD_fit[:,0] = times
        self.MSD_fit[:,1] = MSDs
# =============================================================================
#         print("Full length MSD:", self.length_MSD)
#         print("Fitted length MSD (60 %):", len(MSDs))
# =============================================================================
        try:
            [self.r, self.tau], self.cov = curve_fit(self.function_full_MSD, times, MSDs, p0 = (self.r, self.tau), method = "lm")
        except:
            self.analyse_successful = False
        if self.analyse_successful:
            [chisq, p] = chisquare(MSDs, self.function_full_MSD(times, self.r, self.tau))
            self.chi_MSD_fit = chisq
            self.D_conf = self.r**2.0/(3.0*self.tau)
            self.dr = self.cov[0,0]
            self.dtau = self.cov[1,1]
            self.MSD_fit[:,2] = self.function_full_MSD(times, self.r, self.tau)
            self.MSD_fit[:,3] = self.MSD_fit[:,1] - self.MSD_fit[:,2]
# =============================================================================
#             print(self.MSD_fit)
#             #[chisq, p] = chisquare(MSDs, fit_fn(times))
#         print("Analysis succes?", self.analyse_succesfull)        
# =============================================================================
        
    def plot_full_MSD(self):
        x1, x2 = 0, self.MSD_fit[:,0].max()  # x1 = min, x2 = max
        sp1_y1, sp1_y2 = 0, self.MSD_fit[:,1].max()
        sp2_y1, sp2_y2 = self.MSD_fit[:,3].min(), self.MSD_fit[:,3].max()
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
        sp_1.axis((x1, x2, sp1_y1, sp1_y2))
        sp_2 = plt.subplot2grid((4,4), (3,0), colspan=4, rowspan=1)
        #sp_2 = fig.add_subplot(2, 1, 2) 
        residue_line = np.zeros(len(self.MSD_fit[:,0]))
        sp_2.plot(self.MSD_fit[:,0], residue_line, ":", color = "0.75")
        sp_2.plot(self.MSD_fit[:,0], self.MSD_fit[:,3], "*", color = "0.5", label="residues")
        sp_2.legend()
        sp_2.set_ylabel("Residue")
        sp_2.set_xlabel("Time step [s]")  # Number of data points used in MJD calculation
        sp_2.axis((x1, x2, sp2_y1, sp2_y2))
        plt.show() 
        
    def show_trajectory(self):
        #plt.rc('text', usetex=True)
        #plt.rc('font', family='serif')
        #x1, x2 = self.localizations[:,2].min(), self.localizations[:,2].max()
        #y1, y2 = self.localizations[:,3].min(), self.localizations[:,3].max()
        plt.plot(self.localizations[:,2], self.localizations[:,3], linestyle="--", marker="o")
        plt.title("Trajectory of one particle")
        plt.xlabel(r"$\mu$" + "m in x direction")
        plt.ylabel(r"$\mu$" + "m in y direction")
        #plt.xlabel(r'\textbf{time} (s)')
        plt.show()

    def analyse_particle(self):
        self.calc_trajectory_number()
        self.calc_length_trajectory()
        self.calc_length_MSD()
        self.calc_MSD()
        self.calc_diffusion()
        self.check_immobility()
        if not (self.immobility):
            self.fit_full_MSD()
            self.check_confined()

    def print_particle(self):
        print("Diffusion coefficient:", self.D)
        print("chi2 linear fit:", self.chi_D)
        print("Type immobile:", self.immobility)
        if not self.immobility:
            print("Analyse successful?", self.analyse_successful)
            print("Type confined:", self.confined)
            print("D_conf:", self.D_conf)
            print("r_conf:", self.r)
            print("tau:", self.tau)
            print("tau threshold:", self.tau_threshold)
            print("chi2 rossier fit:", self.chi_MSD_fit)

    def plot_particle(self):
        self.show_trajectory()
        self.plot_diffusion()
        if not self.immobility:
            self.plot_full_MSD()
        elif self.immobility:
            self.plot_full_MSD_immob()
        self.print_particle()

    def save_times_MSDs_60(self):
        out_file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\pySPT_cell_1_MMStack_Pos0\\preAnalysis\\MSD60" + "_" + str(self.molecule_number) + ".txt"
        #out_file_name = directory + "\ " + year + month + day + "_" + base_name + "_trc_format.txt"
        header = "time step [s]\t MSD\t Diffusion coeff: " + str(self.D)
        np.savetxt(out_file_name, 
                   X=self.MSD_60,
                   fmt = ("%.3f","%.14f"),
                   header = header)
    
   
def main():
    trajectory = Trajectory()
    trajectory.analyse_particle()
    trajectory.plot_particle()
        
    
if __name__ == "__main__":
    main()
    