import numpy as np
import math

class Params:
    def __init__(self):
        #Given Parameters
        self.gm1 = 1.02737e5
        self.gm2 = 1.02737e4
        self.bin_sep = 32.6823

        #Derived Parameters
        self.orbital_Omega = np.sqrt((self.gm1+self.gm2)/(self.bin_sep**3)) #1.79925
        self.binary_period = (2*np.pi) / self.orbital_Omega
        self.gamma = 5.0/3.0
        self.mass_ratio = self.gm2/self.gm1
        self.three_one_res = self.bin_sep / ((9*(1+self.mass_ratio))**(1/3))

        #file based properties
        #spherical: self.timesteps_per_filenum = 0.05
        self.timesteps_per_filenum = 0.1
        self.filenums_per_orbit = math.ceil(self.binary_period / self.timesteps_per_filenum)