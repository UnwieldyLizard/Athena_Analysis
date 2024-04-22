#Rust
import rusty_athena
#Code base
from .params import *
from .file_config import*
from .vectors import *
#Logistical
import os
import logging
import warnings
import h5py
from datetime import datetime
#Math
import numpy as np
from scipy.interpolate import griddata
#Matplotlib
from matplotlib import pyplot as plt
from matplotlib import colors as colors
from matplotlib.colors import SymLogNorm, LogNorm
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib import cm as cm
import matplotlib.gridspec as gridspec
from matplotlib import patches
import matplotlib.ticker as mtick
import matplotlib.style as mplstyle
mplstyle.use('fast')
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc #
rc("font", family="DejaVu Sans", weight="normal", size="9") # bryance used Liberation Serif, swapping to DejaVu Sans cuz lib serif not found
rc("axes", grid=True)
rc("grid", linestyle="--")
rc("xtick", direction="in")
rc("ytick", direction="in")
rc("savefig", format="png", bbox="tight")

#building codebase objects
sim = Params()
file = File()
vec = Vector()

DISABLE_VERTICAL_GRAVITY_IN_CYL = True

class Athena_Analysis():
    def __init__(self, filename, grid_type="Spherical"):
        """
        Initializes the file and basic information.
        """
        self.data = h5py.File(filename, 'r')
        self.time = self.data.attrs["Time"]
        self.NumMeshBlocks = self.data.attrs["NumMeshBlocks"]
        #self.DatasetNames = self.data.attrs["DatasetNames"]
        self.VariableNames = self.data.attrs["VariableNames"]
        #self.NumCycles = self.data["NumCycles"]
        #self.Coordinates = self.data.attrs["Coordinates"]
        
        self.gridtype = grid_type

        self.NO_Z_GRAVITY = True

        self.r = None
        self.phi = None
        self.theta = None

        self.x = None
        self.y = None

        self.x3 = None
        self.y3 = None

        self.z = None

        self.dr = None
        self.dtheta = None
        self.dphi = None
        self.dz = None

        self.rccfa = None

        self.cell_vol = None

        self.midplane_plotter_grid = None
        self.angular_cmap = None

        self.wd_grav_pot = None
        self.companion_grav_pot = None
        self.accel_pot = None

        self.total_mass_flux = None # indicator for all mass fluxes

        self.r_vec = None
        self.r_hat = None
        self.v_vec = None
             
    def _axes(self):
        """
        Initializes the base axes.
        """
        if self.gridtype == "Spherical":
            self.r_primitive = np.array(self.data["x1v"], dtype="float64")
            self.theta_primitive = np.array(self.data["x2v"], dtype="float64")
            self.phi_primitive = np.array(self.data["x3v"], dtype="float64")
            self.rf_primitive = np.array(self.data["x1f"], dtype="float64")
            self.thetaf_primitive = np.array(self.data["x2f"], dtype="float64")
            self.phif_primitive = np.array(self.data["x3f"], dtype="float64")

            self.r_len = len(self.r_primitive[0])
            self.theta_len = len(self.theta_primitive[0])
            self.phi_len = len(self.phi_primitive[0])
            self.array_size = (self.NumMeshBlocks, self.phi_len, self.theta_len, self.r_len)
            self.vector_array_size = (3, self.NumMeshBlocks, self.phi_len, self.theta_len, self.r_len)

            self.possible_r = np.sort(list(set((self.r_primitive).flatten())))
            self.possible_theta = np.sort(list(set((self.theta_primitive).flatten())))
            self.possible_phi = np.sort(list(set((self.phi_primitive).flatten())))
            self.possible_rf = np.sort(list(set((self.rf_primitive).flatten())))
            self.possible_thetaf = np.sort(list(set((self.thetaf_primitive).flatten())))
            self.possible_phif = np.sort(list(set((self.phif_primitive).flatten())))

            self.possible_vars = {"r": self.possible_r, "phi": self.possible_phi, "theta": self.possible_theta}

            if len(self.possible_r) == (len(self.possible_rf - 1)):
                self.possible_dr_primitive = self.possible_rf[1:] - self.possible_rf[:-1]
            else:
                self.possible_dr_primitive = np.zeros(len(self.possible_r))
                self.possible_dr_primitive[0] = self.possible_r[0] - 1
                self.possible_dr_primitive[1:] = self.possible_r[1:] - self.possible_r[:-1]
            if len(self.possible_theta) == (len(self.possible_thetaf - 1)):
                self.possible_dtheta_primitive = self.possible_thetaf[1:] - self.possible_thetaf[:-1]
            else:
                self.possible_dtheta_primitive = np.zeros(len(self.possible_theta))
                self.possible_dtheta_primitive[0] = self.possible_theta[0] - 1
                self.possible_dtheta_primitive[1:] = self.possible_theta[1:] - self.possible_theta[:-1]
            if len(self.possible_phi) == (len(self.possible_phif) - 1):
                self.possible_dphi_primitive = self.possible_phif[1:] - self.possible_phif[:-1]
            else:
                self.possible_dphi_primitive = np.zeros(len(self.possible_phi))
                self.possible_dphi_primitive[0] = self.possible_phi[0] - 1
                self.possible_dphi_primitive[1:] = self.possible_phi[1:] - self.possible_phi[:-1]

        if self.gridtype == "Cylindrical":
            self.r_primitive = np.array(self.data["x1v"], dtype="float64")
            self.phi_primitive = np.array(self.data["x2v"], dtype="float64")
            self.z_primitive = np.array(self.data["x3v"], dtype="float64")
            self.rf_primitive = np.array(self.data["x1f"], dtype="float64")
            self.phif_primitive = np.array(self.data["x2f"], dtype="float64")
            self.zf_primitive = np.array(self.data["x3f"], dtype="float64")

            self.r_len = len(self.r_primitive[0])
            self.phi_len = len(self.phi_primitive[0])
            self.z_len = len(self.z_primitive[0])
            self.array_size = (self.NumMeshBlocks, self.z_len, self.phi_len, self.r_len)
            self.vector_array_size = (3, self.NumMeshBlocks, self.z_len, self.phi_len, self.r_len)
            self.tensor_array_size = (3, 3, self.NumMeshBlocks, self.z_len, self.phi_len, self.r_len)

            self.possible_r = np.sort(list(set((self.r_primitive).flatten())))
            self.possible_phi = np.sort(list(set((self.phi_primitive).flatten())))
            self.possible_z = np.sort(list(set((self.z_primitive).flatten())))
            self.possible_rf = np.sort(list(set((self.rf_primitive).flatten())))
            self.possible_phif = np.sort(list(set((self.phif_primitive).flatten())))
            self.possible_zf = np.sort(list(set((self.zf_primitive).flatten())))

            self.possible_vars = {"r": self.possible_r, "phi": self.possible_phi, "z": self.possible_z}

            if len(self.possible_r) == (len(self.possible_rf - 1)):
                self.possible_dr_primitive = self.possible_rf[1:] - self.possible_rf[:-1]
            else:
                self.possible_dr_primitive = np.zeros(len(self.possible_r))
                self.possible_dr_primitive[0] = 2*(self.possible_r[0] - self.possible_rf[0])
                self.possible_dr_primitive[1:] = self.possible_r[1:] - self.possible_r[:-1]
            if len(self.possible_phi) == (len(self.possible_phif - 1)):
                self.possible_dphi_primitive = self.possible_phif[1:] - self.possible_phif[:-1]
            else:
                self.possible_dphi_primitive = np.zeros(len(self.possible_phi))
                self.possible_dphi_primitive[0] = 2*(self.possible_phi[0] - self.possible_phif[0])
                self.possible_dphi_primitive[1:] = self.possible_phi[1:] - self.possible_phi[:-1]
            if len(self.possible_z) == (len(self.possible_zf) - 1):
                self.possible_dz_primitive = self.possible_zf[1:] - self.possible_zf[:-1]
            else:
                self.possible_dz_primitive = np.zeros(len(self.possible_z))
                self.possible_dz_primitive[0] = 2*(self.possible_z[0] - self.possible_zf[0])
                self.possible_dz_primitive[1:] = self.possible_z[1:] - self.possible_z[:-1]
            
    def get_primaries(self, get_rho=False, get_press=False, get_vel_r=False, get_vel_theta=False, get_vel_phi=False, get_vel_z=False):
        """
        Generates arrays of primary data.
        Athena data arrays are always 4 dimensional with axis 0 interating through the mesh blocks, axis 1 through phi, axis 2 through theta, and axis 3 through r.

        Parameters
        ----------
        get_rho : bool, default False
            Determines if rho array is generated.
        get_press : bool, default False
            Determines if press array is generated.
        get_vel1 : bool, default False
            Determines if vel1 array is generated.
        get_vel2 : bool, default False
            Determines if vel2 array is generated.
        get_vel3 : bool, default False
            Determines if vel3 array is generated.
        """
        if get_rho == True:
            self.rho = np.array(self.data["prim"][0])
        if get_press == True:
            self.press = np.array(self.data["prim"][1])
        if self.gridtype == "Spherical":
            if get_vel_r == True:
                self.vel_r = np.array(self.data["prim"][2])
            if get_vel_theta == True:
                self.vel_theta = np.array(self.data["prim"][3])
            if get_vel_phi == True:
                self.vel_phi = np.array(self.data["prim"][4])
            if get_vel_z == True:
                warnings.warn("Inappropirate coordinates: Calling for z in a spherical geometry")
        if self.gridtype == "Cylindrical":
            if get_vel_r == True:
                self.vel_r = np.array(self.data["prim"][2])
            if get_vel_phi == True:
                self.vel_phi = np.array(self.data["prim"][3])
            if get_vel_z == True:
                self.vel_z = np.array(self.data["prim"][4])
            if get_vel_theta == True:
                warnings.warn("Inappropirate coordinates: Calling for theta in a cylindrical geometry. Perhaps you meant phi.")

    def get_Bfields(self):
        """
        Gets the components of the B field
        """
        if self.gridtype == "Spherical":
            self.B_r = self.data["B"][0] #Bcc1
            self.B_theta = self.data["B"][1] #Bcc2
            self.B_phi = self.data["B"][2] #Bcc3
        if self.gridtype == "Cylindrical":
            self.B_r = self.data["B"][0] #Bcc1
            self.B_phi = self.data["B"][1] #Bcc2
            self.B_z = self.data["B"][2] #Bcc3

    def spherical_grid(self, get_r=False, get_theta=False, get_phi=False, get_dr=False, get_dtheta=False, get_dphi=False):
        """
        DEPRECIATED PLEASE SWITCH TO NATIVE GRID
        Generates full grid of simulation space in spherical coordinates.

        Parameters
        ----------
        get_r : bool, default False
            Determines if r axis is generated.
        get_theta : bool, default False 
            Determines if theta axis is generated.
        get_phi: bool, default False
            Determines if phi axis is generated.
        """
        self._axes()
        if get_r == True and self.r is None:
            self.r = np.array([[[self.r_primitive[k] for i in range(self.theta_len)] for j in range(self.phi_len)]for k in range(self.NumMeshBlocks)])
        if get_theta == True and self.theta is None:
            self.theta = np.array([[[self.theta_primitive[k] for i in range(self.r_len)] for j in range(self.phi_len)]for k in range(self.NumMeshBlocks)])
            self.theta = np.swapaxes(self.theta,2,3)
        if get_phi == True and self.phi is None: 
            self.phi = np.array([[[self.phi_primitive[k] for i in range(self.theta_len)] for j in range(self.r_len)]for k in range(self.NumMeshBlocks)])
            self.phi = np.swapaxes(self.phi,1,3)

        if get_dr == True and self.dr is None:
            self.dr = np.zeros((self.NumMeshBlocks, self.r_len))
            for n in range(self.NumMeshBlocks):
                self.dr[n] = self.rf_primitive[n, 1:] - self.rf_primitive[n, :-1]
            self.dr = np.array([[[self.dr[k] for i in range(self.theta_len)] for j in range(self.phi_len)]for k in range(self.NumMeshBlocks)])   

        if get_dtheta == True and self.dtheta is None:
            self.dtheta = np.zeros((self.NumMeshBlocks, self.theta_len))
            for n in range(self.NumMeshBlocks):
                self.dtheta[n] = self.thetaf_primitive[n, 1:] - self.thetaf_primitive[n, :-1]
            self.dtheta = np.array([[[self.dtheta[k] for i in range(self.r_len)] for j in range(self.phi_len)]for k in range(self.NumMeshBlocks)])
            self.dtheta = np.swapaxes(self.dtheta,2,3)

        if get_dphi == True and self.dphi is None:
            self.dphi = np.zeros((self.NumMeshBlocks, self.phi_len))
            for n in range(self.NumMeshBlocks):
                self.dphi[n] = self.phif_primitive[n, 1:] - self.phif_primitive[n, :-1]
            self.dphi = np.array([[[self.dphi[k] for i in range(self.theta_len)] for j in range(self.r_len)]for k in range(self.NumMeshBlocks)])
            self.dphi = np.swapaxes(self.dphi,1,3)

    def native_grid(self, get_r=False, get_theta=False, get_phi=False, get_z=False, get_dr=False, get_dtheta=False, get_dphi=False, get_dz=False):
            """
            Generates full grid of simulation space in it's native coordinates.

            Parameters
            ----------
            get_r : bool, default False
                Determines if r axis is generated.
            get_theta : bool, default False 
                Determines if theta axis is generated.
            get_phi: bool, default False
                Determines if phi axis is generated.
            get_z : bool, default False
                Determines if z axis is generated.
            get_dr : bool, default False
                Determines if r differential is generated.
            get_dtheta : bool, default False 
                Determines if theta differential is generated.
            get_dphi: bool, default False
                Determines if phi differential is generated.
            get_dz : bool, default False
                Determines if z differential is generated.
            """
            if self.gridtype == "Spherical":
                self._axes()
                if get_r == True and self.r is None:
                    self.r = np.array([[[self.r_primitive[k] for i in range(self.theta_len)] for j in range(self.phi_len)]for k in range(self.NumMeshBlocks)])
                if get_theta == True and self.theta is None:
                    self.theta = np.array([[[self.theta_primitive[k] for i in range(self.r_len)] for j in range(self.phi_len)]for k in range(self.NumMeshBlocks)])
                    self.theta = np.swapaxes(self.theta,2,3)
                if get_phi == True and self.phi is None: 
                    self.phi = np.array([[[self.phi_primitive[k] for i in range(self.theta_len)] for j in range(self.r_len)]for k in range(self.NumMeshBlocks)])
                    self.phi = np.swapaxes(self.phi,1,3)

                if get_dr == True and self.dr is None:
                    self.dr = np.zeros((self.NumMeshBlocks, self.r_len))
                    for n in range(self.NumMeshBlocks):
                        self.dr[n] = self.rf_primitive[n, 1:] - self.rf_primitive[n, :-1]
                    self.dr = np.array([[[self.dr[k] for i in range(self.theta_len)] for j in range(self.phi_len)]for k in range(self.NumMeshBlocks)])   

                if get_dtheta == True and self.dtheta is None:
                    self.dtheta = np.zeros((self.NumMeshBlocks, self.theta_len))
                    for n in range(self.NumMeshBlocks):
                        self.dtheta[n] = self.thetaf_primitive[n, 1:] - self.thetaf_primitive[n, :-1]
                    self.dtheta = np.array([[[self.dtheta[k] for i in range(self.r_len)] for j in range(self.phi_len)]for k in range(self.NumMeshBlocks)])
                    self.dtheta = np.swapaxes(self.dtheta,2,3)

                if get_dphi == True and self.dphi is None:
                    self.dphi = np.zeros((self.NumMeshBlocks, self.phi_len))
                    for n in range(self.NumMeshBlocks):
                        self.dphi[n] = self.phif_primitive[n, 1:] - self.phif_primitive[n, :-1]
                    self.dphi = np.array([[[self.dphi[k] for i in range(self.theta_len)] for j in range(self.r_len)]for k in range(self.NumMeshBlocks)])
                    self.dphi = np.swapaxes(self.dphi,1,3)
                
                if get_z == True or get_dz == True:
                    warnings.warn("Inappropirate coordinates: Calling for z in a spherical geometry")
            if self.gridtype == "Cylindrical":
                self._axes()
                if get_r == True and self.r is None:
                    self.r = np.array([[[self.r_primitive[k] for i in range(self.phi_len)] for j in range(self.z_len)]for k in range(self.NumMeshBlocks)])
                if get_phi == True and self.phi is None:
                    self.phi = np.array([[[self.phi_primitive[k] for i in range(self.r_len)] for j in range(self.z_len)]for k in range(self.NumMeshBlocks)])
                    self.phi = np.swapaxes(self.phi,2,3)
                if get_z == True and self.z is None: 
                    self.z = np.array([[[self.z_primitive[k] for i in range(self.phi_len)] for j in range(self.r_len)]for k in range(self.NumMeshBlocks)])
                    self.z = np.swapaxes(self.z,1,3)

                if get_dr == True and self.dr is None:
                    self.dr = np.zeros((self.NumMeshBlocks, self.r_len))
                    for n in range(self.NumMeshBlocks):
                        self.dr[n] = self.rf_primitive[n, 1:] - self.rf_primitive[n, :-1]
                    self.dr = np.array([[[self.dr[k] for i in range(self.phi_len)] for j in range(self.z_len)]for k in range(self.NumMeshBlocks)])

                if get_dphi == True and self.dphi is None:
                    self.dphi = np.zeros((self.NumMeshBlocks, self.phi_len))
                    for n in range(self.NumMeshBlocks):
                        self.dphi[n] = self.phif_primitive[n, 1:] - self.phif_primitive[n, :-1]
                    self.dphi = np.array([[[self.dphi[k] for i in range(self.r_len)] for j in range(self.z_len)]for k in range(self.NumMeshBlocks)])
                    self.dphi = np.swapaxes(self.dphi,2,3)

                if get_dz == True and self.dz is None:
                    self.dz = np.zeros((self.NumMeshBlocks, self.z_len))
                    for n in range(self.NumMeshBlocks):
                        self.dz[n] = self.zf_primitive[n, 1:] - self.zf_primitive[n, :-1]
                    self.dz = np.array([[[self.dz[k] for i in range(self.phi_len)] for j in range(self.r_len)]for k in range(self.NumMeshBlocks)])
                    self.dz = np.swapaxes(self.dz,1,3)

                if get_theta == True or get_dtheta == True:
                    warnings.warn("Inappropirate coordinates: Calling for theta in a cylindrical geometry. Perhaps you meant phi.")

    def get_cell_vol(self):
        """
        Gets the volumes of all cells. Calculaton does using integral over one cell in native coordinates.
        """
        if self.cell_vol is None and self.gridtype == "Spherical":
            self._axes()
            self.cell_vol = np.zeros((self.NumMeshBlocks, self.phi_len, self.theta_len, self.r_len))
            for i in range(self.NumMeshBlocks):
                rpart = (self.rf_primitive[i][1:]**3 - self.rf_primitive[i][:-1]**3)/3.
                tpart = -(np.cos(self.thetaf_primitive[i][1:]) - np.cos(self.thetaf_primitive[i][:-1]))
                ppart = (self.phif_primitive[i][1:] - self.phif_primitive[i][:-1])
                self.cell_vol[i] = np.einsum("k,j,i->kji", ppart, tpart, rpart)

        if self.cell_vol is None and self.gridtype == "Cylindrical":
            self._axes()
            self.cell_vol = np.zeros((self.NumMeshBlocks, self.phi_len, self.theta_len, self.r_len))
            for i in range(self.NumMeshBlocks):
                rpart = (self.rf_primitive[i][1:]**2 - self.rf_primitive[i][:-1]**2)/2.
                ppart = (self.phif_primitive[i][1:] - self.phif_primitive[i][:-1])
                zpart = (self.zf_primitive[i][1:] - self.zf_primitive[i][:-1])
                self.cell_vol[i] = np.einsum("k,j,i->kji", zpart, ppart, rpart)

    def get_face_areas(self, get_rcc_face_areas = False, get_r_face_areas = False):
        """
        Gets r cell centered face areas.
        """
        self._axes()
        if get_rcc_face_areas == True and self.rccfa is None and self.gridtype == "Spherical":
            self.rcc_face_area = np.zeros(self.array_size)
            for i in range(self.NumMeshBlocks):
                rpart = self.r_primitive[i] ** 2
                tpart = -(np.cos(self.thetaf_primitive[i][1:]) - np.cos(self.thetaf_primitive[i][:-1]))
                ppart = (self.phif_primitive[i][1:] - self.phif_primitive[i][:-1])
                self.rcc_face_area[i] = np.einsum("k,j,i->kji", ppart, tpart, rpart)

        if get_rcc_face_areas == True and self.rccfa is None and self.gridtype == "Cylindrical":
            self.rcc_face_area = np.zeros(self.array_size)
            for i in range(self.NumMeshBlocks):
                rpart = self.r[i] # self.r_primitive[i]
                ppart = self.dphi[i] #(self.phif_primitive[i][1:] - self.phif_primitive[i][:-1])
                zpart = self.dz[i] #(self.zf_primitive[i][1:] - self.zf_primitive[i][:-1])
                self.rcc_face_area[i] = rpart * ppart * zpart # np.einsum("k,j,i->kji", zpart, ppart, rpart)

        self._axes()
        if get_r_face_areas == True and self.gridtype == "Spherical":
            self.r_face_area = np.zeros(self.array_size + np.array([0,0,0,1])) # adding 1 to radial cuz one more faces than cells
            for i in range(self.NumMeshBlocks):
                rpart = self.rf_primitive[i] ** 2
                tpart = -(np.cos(self.thetaf_primitive[i][1:]) - np.cos(self.thetaf_primitive[i][:-1]))
                ppart = (self.phif_primitive[i][1:] - self.phif_primitive[i][:-1])
                self.r_face_area[i] = np.einsum("k,j,i->kji", ppart, tpart, rpart)

        if get_r_face_areas == True and self.gridtype == "Cylindrical":
            self.r_face_area = np.zeros(self.array_size + np.array([0,0,0,1])) # adding 1 to radial cuz one more faces than cells
            for i in range(self.NumMeshBlocks):
                rpart = self.rf_primitive[i]
                ppart = (self.phif_primitive[i][1:] - self.phif_primitive[i][:-1])
                zpart = (self.zf_primitive[i][1:] - self.zf_primitive[i][:-1])
                self.r_face_area[i] = np.einsum("k,j,i->kji", zpart, ppart, rpart)

    def cart_grid(self, projection):
        """
        Generates full grid of simulation space in cartesian coordinates. The cooresponding native grid will also be generated as an intermediate step if it does not already exist.

        Parameters
        ----------
        projection : str
            If '2D' only generates 2D cartesian grid of x and y in the midplane. If '3D' generates 3D cartesian grid of entire simulation space.
        """
        if self.gridtype == "Spherical":
        #if (projection == '2D' or projection == '2d') and (self.x is None or self.y is None):
            self.spherical_grid(get_r=True, get_phi=True)
            if self.x is None:
                self.x = self.r * np.cos(self.phi)
            if self.y is None:
                self.y = self.r * np.sin(self.phi)
            if (projection == '3D' or projection == '3d') and (self.x is None or self.y is None or self.z is None):
                self.spherical_grid(get_r=True, get_theta=True, get_phi=True)
                R = self.r * np.sin(self.theta)
                if self.x3 is None:
                    self.x3 = R * np.cos(self.phi)
                if self.y3 is None:
                    self.y3 = R * np.sin(self.phi)
                if self.z is None:
                    self.z = self.r * np.cos(self.theta)

        if self.gridtype == "Cylindrical":
            self.native_grid(get_r=True, get_phi=True)
            if self.x is None:
                self.x = self.r * np.cos(self.phi)
            if self.y is None:
                self.y = self.r * np.sin(self.phi)
            if projection == "3D" or projection == '3d':
                self.native_grid(get_z=True)

        if self.gridtype == "Cartesian":
            warning.warn("Just why, the native grid is already cartesian lol")
            self.native_grid(get_x=True, get_y=True, get_z=True)

    def _slice_filter(self, idx):
        """
        Filter function for coordinate slices
        """
        if self.slice_axis == 'r':
            if any(self.slice_lower_bound < ele < self.slice_upper_bound for ele in self.r_primitive[idx]):
                return True
        if self.slice_axis == 'theta':
            if any(self.slice_lower_bound < ele < self.slice_upper_bound for ele in self.theta_primitive[idx]):
                return True
        if self.slice_axis == 'phi':
            if any(self.slice_lower_bound < ele < self.slice_upper_bound for ele in self.phi_primitive[idx]):
                return True
        if self.slice_axis == 'z':
            if any(self.slice_lower_bound < ele < self.slice_upper_bound for ele in self.z_primitive[idx]):
                return True    

    def get_slice_idx(self, slice_axis, slice_center, slice_width):
        """
        Generates array of idices within a coordinate slice which can then be used to filter data arrays for only those within this slice.

        Parameters
        ----------
        slice_axis : str
            Accepts 'r', 'theta', 'phi' , 'x', 'y', or 'z' then performs the slice along that axis.
        slice_center : float
            Center value of the slice.
        slice_width : float
            Width of the slice.
        """
        self._axes()
        if self.gridtype == "Spherical":
            idx = np.arange(0, self.NumMeshBlocks)
            if slice_axis == 'z':
                self.cart_grid(projection='3D')
            if slice_axis == 'r' or slice_axis == 'theta' or slice_axis == 'phi':
                self.slice_axis = slice_axis
                self.slice_lower_bound = slice_center - (slice_width/2)
                self.slice_upper_bound = slice_center + (slice_width/2)
                slice_idx = np.array(list(filter(self._slice_filter, idx)))
                #print(slice_axis, "slice (no y)", self.slice_idx)
        if self.gridtype == "Cylindrical":
            idx = np.arange(0, self.NumMeshBlocks)
            if slice_axis == 'z':
                self.cart_grid(projection='3D')
            if slice_axis == 'r' or slice_axis == 'phi': # or slice_axis == 'z'
                #this don't work for z and idk why
                self.slice_axis = slice_axis
                self.slice_lower_bound = slice_center - (slice_width/2)
                self.slice_upper_bound = slice_center + (slice_width/2)
                slice_idx = np.array(list(filter(self._slice_filter, idx)))
        if slice_axis == 'x':
            self.cart_grid(projection='2D')
            slice_idx = [abs(self.x - slice_center) < slice_width]
        if slice_axis == 'y':
            self.cart_grid(projection='2D')
            slice_idx = [abs(self.y - slice_center) < slice_width]
        if slice_axis == 'z':
            #this works and idk why
            self.cart_grid(projection='2D')
            slice_idx = [abs(self.z - slice_center) < slice_width]
        if len(slice_idx) == 0:
            warnings.warn("The slice is empty. This will raise errors if you try to index with it.")
        return slice_idx

    def get_rot_slice_idx(self, rotation, slice_center, slice_width):
        self._axes()
        if self.gridtype == "Cylindrical":
            self.cart_grid(projection='3D')
            idx = np.arange(0, self.NumMeshBlocks)
            complex_position = (self.x + self.y*1j) * np.exp(-1j*rotation)
            slice_x = complex_position.real
            slice_y = complex_position.imag
            slice_idx = [abs(slice_y - slice_center) < slice_width]
        if len(slice_idx) == 0:
            warnings.warn("The slice is empty. This will raise errors if you try to index with it.")
        return slice_idx

    def _midplane_plotter_grid(self, slice_xmax, slice_ymax):
        [slice_x, slice_y] = (np.mgrid[-slice_xmax:slice_xmax:512j,-slice_ymax:slice_ymax:512j])
        return [slice_x, slice_y]

    def _build_angular_cmap(self):
        if self.angular_cmap is None:
            newcolors = np.zeros((256, 4))
            for n in range(256):
                #set red
                if n < 128:
                    #set red
                    newcolors[n][0] = (np.sin(np.linspace(0, np.pi, 256)) ** 1)[n+128]
                elif n >= 128 and n < 192:
                    newcolors[n][0] = 0
                elif n >= 192:
                    newcolors[n][0] = (np.sin(np.linspace(0, np.pi, 128)) ** 1)[n-192]
                #set green
                if n < 64:
                    newcolors[n][1] = (np.sin(np.linspace(0, np.pi, 128)) ** 1)[n+64]
                elif n >= 64 and n < 128:
                    newcolors[n][1] = 0
                elif n >= 128:
                    newcolors[n][1] = (np.sin(np.linspace(0, np.pi, 256)) ** 1)[n-128]
                #set blue
                if n < 256:
                    newcolors[n][2] = (np.sin(np.linspace(0, np.pi, 256)) ** 2)[n]
                #elif n >= 192:
                    #newcolors[n][2] = 0
                #set lightness
                newcolors[n][3] = 1
            self.angular_cmap = colors.ListedColormap(newcolors, "angular")

    def midplane_colorplot(self, q, ax, slicetype='z', log=True, vbound=[1e-5, 1e2], angular=False, plot_COM=False, cmap="viridis", norm=None, rotation = 0):
        """
        Makes a color plot of a value over the midplane

        Parameters
        ----------
        q : list
            The athena data being plotted
        ax : str
            The name of the matplot ax you want to use for this plot
        slicetype : str, default theta
            Determines how the midplane slice is calculated, if 'theta' it will take eveything within a certain theta of theta = pi/2, if 'z' it will take everything within a certain z of z=0.
        log : bool, default True
            Determine if Log scaled axes are used
        vbound : list of int, default [1e-5, 10]
            Determins the bounds of the y axis.
        rotation : float, default 0 
            The angle in radians the image should be rotated (not it does this by rotating the coordinates the opposite amount)
        """

        if angular == True:
            if log == True:
                warnings.warn("log and angular were both true, you shouldn't log angular quantities, setting log to false")
                log=False

        #z slice
        if slicetype == 'z':
            self.cart_grid(projection='3D')
            slice_idx = self.get_slice_idx(slice_axis='z', slice_center=0, slice_width=0.6) #was 0.6

            #print('z', q.shape)
            q = q[tuple(slice_idx)]
            #print('zs', q.shape)
            x = self.x[tuple(slice_idx)]
            y = self.y[tuple(slice_idx)]
            complex_position = (x + y*1j) * np.exp(-1j*rotation)
            x = complex_position.real
            y = complex_position.imag
            z = self.z[tuple(slice_idx)]
            ax.set_xlabel("x")
            ax.set_ylabel("y")
        #z slice end

        #y slice profile
        if slicetype == 'y':
            self.cart_grid(projection='3D')
            slice_idx = self.get_slice_idx(slice_axis='y', slice_center=0, slice_width=0.6)
        
            #print("y", q.shape)
            q = q[tuple(slice_idx)]
            #print("ys", q.shape)
            x = self.x[tuple(slice_idx)]
            y = self.y[tuple(slice_idx)]
            if self.gridtype == "Cylindrical":
                z = self.z[tuple(slice_idx)]
            else:
                z = self.z[tuple(slice_idx)]

            ax.set_xlabel("x")
            ax.set_ylabel("z")
        #y slice end

        #vert slice profile
        if slicetype == 'vert':
            self.cart_grid(projection='3D')
            slice_idx = self.get_rot_slice_idx(rotation=rotation, slice_center=0, slice_width=1)
        
            #print("y", q.shape)
            q = q[tuple(slice_idx)]
            #print("ys", q.shape)
            x = self.x[tuple(slice_idx)]
            y = self.y[tuple(slice_idx)]
            complex_position = (x + y*1j) * np.exp(-1j*rotation)
            x = complex_position.real
            y = complex_position.imag
            if self.gridtype == "Cylindrical":
                z = self.z[tuple(slice_idx)]
            else:
                z = self.z[tuple(slice_idx)]

            ax.set_xlabel("x")
            ax.set_ylabel("z")
        #vert slice end

        #theta slice
        if slicetype == 'theta':
            self.cart_grid(projection='2D')
            slice_idx = self.get_slice_idx(slice_axis='theta', slice_center=np.pi/2, slice_width=0.005*np.pi) #for CVThin2 0.0014*pi was used

            x = self.x[slice_idx]
            y = self.y[slice_idx]
            #print("t", q.shape)
            q = q[slice_idx]
            #print("ts", q.shape)

            ax.set_xlabel("x")
            ax.set_ylabel("y")
        #theta slice end

        cone_theta = np.pi/2

        #theta_idx = np.argmin(np.abs(theta[0, :, 0] - cone_theta))

        # cone
        self.native_grid(get_r=True, get_z=True)
        if self.gridtype == 'Cylindrical' and (slicetype == 'y' or slicetype == 'vert'):
            slice_xmax = self.possible_rf[-1]
            slice_ymax = self.possible_zf[-1]
        else:
            slice_xmax = self.possible_rf[-1]
            slice_ymax = self.possible_rf[-1]

        [slice_x, slice_y] = self._midplane_plotter_grid(slice_xmax, slice_ymax)
            
        if slicetype == 'z' or slicetype == 'theta':
            point = np.array([x.reshape((-1)), y.reshape((-1))], dtype=np.float32).T
        if slicetype == 'y':
            point = np.array([x.reshape((-1)), z.reshape((-1))], dtype=np.float32).T
        if slicetype == "vert":
            point = np.array([x.reshape((-1)), z.reshape((-1))], dtype=np.float32).T
        q_midplane = griddata(point, q.reshape((-1)), (slice_x, slice_y), method="nearest")    

        if vbound == [] or vbound == 0:
            vmin = np.min(q)
            vmax = np.max(q)
        else:
            if isinstance(vbound, list):
                [vmin, vmax] = vbound
            else:
                vmax = np.abs(vbound)
                vmin = -vmax

        clip=True

        if norm is None:
            if log == True and vmax > 0 and vmin < 0:
                norm = SymLogNorm(linthresh=0.01, linscale=1, vmin=vmin, vmax=vmax)
            elif log == True: 
                norm = colors.LogNorm(vmin, vmax, clip)
            else:
                norm = colors.Normalize(vmin, vmax, clip)

        # circle for clipping beyond outer boundary
        if self.gridtype == 'Cylindrical' and (slicetype == 'y' or slicetype == 'vert'):
            outer_bound = patches.Rectangle((-slice_xmax, -slice_ymax), 2*slice_xmax, 2*slice_ymax, facecolor='none')
        else:
            outer_bound = patches.Circle((0, 0), radius = float(self.data.attrs["RootGridX1"][1]), facecolor='none')
        ax.add_patch(outer_bound)

        #special colormaps
        if log == True and vmin < 0 and vmax > 0:
            cmap="PRGn"
        if angular == True:
            self._build_angular_cmap()
            cmap=self.angular_cmap
            
        ax.grid(False)
        if angular == True:
            q_midplane = q_midplane / np.pi

        im = ax.pcolormesh(slice_x.T, slice_y.T, q_midplane.T, cmap=cmap, norm=norm, clip_path=outer_bound, clip_on=True, shading="auto") # readd after norm: clip_path=outer_bound
        ax.set_title(r"%s" % (""))
        if self.gridtype == 'Cylindrical' and (slicetype == 'y' or slicetype == 'vert'):
            ax.set_aspect('auto')
        else:
            ax.set_aspect(1)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        ax.grid(False)
        cbar = plt.colorbar(im, cax=cax)
        cbar.ax.tick_params()
        if angular == True:
            cbar.ax.set_yticklabels(['0',r'$\frac{\pi}{2}$',r'$\pi$',r'$\frac{3\pi}{2}$',r'2$\pi$'])

        #add the white dwarf
        wd_radius=1
        if self.gridtype == 'Cylindrical' and (slicetype == 'y' or slicetype == 'vert'):
            wd = plt.Rectangle((-wd_radius, -slice_ymax), wd_radius*2, slice_ymax*2, color="k")
        else:
            wd = plt.Circle((0, 0), wd_radius*np.sin(cone_theta), color="k")
        ax.add_artist(wd)

        #add the resonant radius
        #res_radius=15.2
        res_radius = sim.three_one_res
        if (slicetype == 'y' or slicetype == 'vert'):
            res_circ = plt.Rectangle((-res_radius, -slice_ymax), res_radius*2, slice_ymax*2, color="w", fill=False, linewidth=0.5)
        else:
            res_circ = plt.Circle((0, 0), res_radius*np.sin(cone_theta), color="w", fill=False, linewidth=0.5)
        ax.add_artist(res_circ)

        #add center of mass marker
        if plot_COM == True:
            self.cart_grid(projection = '2D')
            total_mass = self.integrate(self.rho, 'all')
            cm_x = self.integrate(self.rho * self.x, 'all')
            cm_y = self.integrate(self.rho * self.y, 'all')
            [cm_x, cm_y] = [cm_x /total_mass, cm_y / total_mass]
            cm = plt.Circle((cm_x, cm_y), 0.5, color=[1, 0, 0.7071])
            ax.add_artist(cm)

        ax.set_xlim([-slice_xmax, slice_xmax])
        ax.set_ylim([-slice_ymax, slice_ymax])
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2g'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))

    def midplane_vectorplot(self, q_0, ax, slicetype='z', log=True, vbound=[1e-5, 1e2], rotation = 0, cmap="plasma"):
        """
        Makes a color plot of a value over the midplane

        Parameters
        ----------
        q : list
            The athena data being plotted
        ax : str
            The name of the matplot ax you want to use for this plot
        slicetype : str, default theta
            Determines how the midplane slice is calculated, if 'theta' it will take eveything within a certain theta of theta = pi/2, if 'z' it will take everything within a certain z of z=0.
        log : bool, default True
            Determine if Log scaled axes are used
        vbound : list of int, default [1e-5, 10]
            Determins the bounds of the y axis.
        rotation : float, default 0 
            The angle in radians the image should be rotated (not it does this by rotating the coordinates the opposite amount)
        """
        q_0 = np.array(q_0)
        q = [[],[]]


        if (vbound is not None) and (vbound[0] < 0 and log == True):
            log = False
            warnings.warn("disabled log scale since your bounds extend below zero")
        
        #z slice
        if slicetype == 'z':
            self.cart_grid(projection='3D')
            slice_idx = self.get_slice_idx(slice_axis='z', slice_center=0, slice_width=0.6) #was 0.6

            #print('z', q.shape)
            for i in range(q_0.shape[0]):
                q[i] = q_0[i][tuple(slice_idx)]
            #print('zs', q.shape)
            x = self.x[tuple(slice_idx)]
            y = self.y[tuple(slice_idx)]
            complex_position = (x + y*1j) * np.exp(-1j*rotation)
            x = complex_position.real
            y = complex_position.imag
            z = self.z[tuple(slice_idx)]
            ax.set_xlabel("x")
            ax.set_ylabel("y")
        #z slice end

        #y slice profile
        if slicetype == 'y':
            self.cart_grid(projection='3D')
            slice_idx = self.get_slice_idx(slice_axis='y', slice_center=0, slice_width=0.6)
        
            #print("y", q.shape)
            for i in range(range(q_0.shape[0])):
                q[i] = q_0[i][tuple(slice_idx)]
            #print("ys", q.shape)
            x = self.x[tuple(slice_idx)]
            y = self.y[tuple(slice_idx)]
            if self.gridtype == "Cylindrical":
                z = self.z[tuple(slice_idx)]
            else:
                z = self.z[tuple(slice_idx)]

            ax.set_xlabel("x")
            ax.set_ylabel("z")
        #y slice end

        #vert slice profile
        if slicetype == 'vert':
            self.cart_grid(projection='3D')
            slice_idx = self.get_rot_slice_idx(rotation=rotation, slice_center=0, slice_width=1)
        
            #print("y", q.shape)
            for i in range(q_0.shape[0]):
                q[i] = q_0[i][tuple(slice_idx)]
            #print("ys", q.shape)
            x = self.x[tuple(slice_idx)]
            y = self.y[tuple(slice_idx)]
            complex_position = (x + y*1j) * np.exp(-1j*rotation)
            x = complex_position.real
            y = complex_position.imag
            if self.gridtype == "Cylindrical":
                z = self.z[tuple(slice_idx)]
            else:
                z = self.z[tuple(slice_idx)]

            ax.set_xlabel("x")
            ax.set_ylabel("z")
        #vert slice end

        #theta slice
        if slicetype == 'theta':
            self.cart_grid(projection='2D')
            slice_idx = self.get_slice_idx(slice_axis='theta', slice_center=np.pi/2, slice_width=0.005*np.pi) #for CVThin2 0.0014*pi was used

            x = self.x[slice_idx]
            y = self.y[slice_idx]
            #print("t", q.shape)
            q = q[:, slice_idx]
            #print("ts", q.shape)

            ax.set_xlabel("x")
            ax.set_ylabel("y")
        #theta slice end

        self.native_grid(get_r=True, get_z=True)
        if self.gridtype == 'Cylindrical' and (slicetype == 'y' or slicetype == 'vert'):
            slice_xmax = self.possible_rf[-1]
            slice_ymax = self.possible_zf[-1]
        else:
            slice_xmax = self.possible_rf[-1]
            slice_ymax = self.possible_rf[-1]

        [slice_x, slice_y] = self._midplane_plotter_grid(slice_xmax, slice_ymax)
            
        if slicetype == 'z' or slicetype == 'theta':
            point = np.array([x.reshape((-1)), y.reshape((-1))], dtype=np.float32).T
        if slicetype == 'y':
            point = np.array([x.reshape((-1)), z.reshape((-1))], dtype=np.float32).T
        if slicetype == "vert":
            point = np.array([x.reshape((-1)), z.reshape((-1))], dtype=np.float32).T
        q_midplane = np.array([griddata(point, q[0].reshape((-1)), (slice_x, slice_y), method="nearest"), griddata(point, q[1].reshape((-1)), (slice_x, slice_y), method="nearest")])    
      
        magnitudes = np.sqrt(q_midplane[0]*q_midplane[0] + q_midplane[1]*q_midplane[1])

        if vbound is None:
            vmin = np.min(magnitudes)
            vmax = np.max(magnitudes)
        else:
            if isinstance(vbound, list):
                [vmin, vmax] = vbound
            else:
                vmax = np.abs(vbound)
                vmin = -vmax

        mask = (magnitudes < (vmin + 0.001*(vmax-vmin)))
        q_midplane[0][mask] = np.nan
        q_midplane[1][mask] = np.nan

        clip=True
        if log == True:
            norm = colors.LogNorm(vmin, vmax, clip)
        else:
            norm = colors.Normalize(vmin, vmax, clip)

        # circle for clipping beyond outer boundary
        if not (self.gridtype == 'Cylindrical' and (slicetype == 'y' or slicetype == 'vert')):
            mask = (slice_x*slice_x+slice_y*slice_y > float(self.data.attrs["RootGridX1"][1])*float(self.data.attrs["RootGridX1"][1]))
            q_midplane[0][mask] = np.nan
            q_midplane[1][mask] = np.nan
            mask = (slice_x*slice_x+slice_y*slice_y < 1)
            q_midplane[0][mask] = np.nan
            q_midplane[1][mask] = np.nan
        
        ax.grid(False)

        im = ax.streamplot(slice_x.T, slice_y.T, q_midplane[0].T, q_midplane[1].T, color=magnitudes.T, cmap=cmap, norm=norm, broken_streamlines=False) # readd after norm: clip_path=outer_bound,  clip_on=True
        #qk = ax.quiverkey(im, 0.9, 0.9, 300, r'$3000$', labelpos='E', coordinates='figure')
        ax.set_title(r"%s" % (""))
        if self.gridtype == 'Cylindrical' and (slicetype == 'y' or slicetype == 'vert'):
            ax.set_aspect('auto')
        else:
            ax.set_aspect(1)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", size="5%", pad=0.05)
        ax.grid(False)
        cbar = plt.colorbar(im.lines, cax=cax, orientation="horizontal")
        cax.xaxis.set_ticks_position("top")
        cbar.ax.tick_params()

        #add the white dwarf
        wd_radius=1
        if self.gridtype == 'Cylindrical' and (slicetype == 'y' or slicetype == 'vert'):
            wd = plt.Rectangle((-wd_radius, -slice_ymax), wd_radius*2, slice_ymax*2, color="k")
        else:
            wd = plt.Circle((0, 0), wd_radius, color="k")
        ax.add_artist(wd)

        #add the resonant radius
        #res_radius=15.2
        res_radius = sim.three_one_res
        if (slicetype == 'y' or slicetype == 'vert'):
            res_circ = plt.Rectangle((-res_radius, -slice_ymax), res_radius*2, slice_ymax*2, color="w", fill=False, linewidth=0.5)
        else:
            res_circ = plt.Circle((0, 0), res_radius, color="w", fill=False, linewidth=0.5)
        ax.add_artist(res_circ)

        ax.set_xlim([-slice_xmax, slice_xmax])
        ax.set_ylim([-slice_ymax, slice_ymax])
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2g'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))


    def _ARCHAIC_azmuthal_integral(self, q):
        """
        Azmuthally integrates an athena data value. Note the returned integrated data will be in an entirely different configuration.
        The returned array will be 2D with axis 0 indexing through theta values and axis 1 indexing through r values.
        To find the r and theta values of a point find the values of the same index in self.r_possible and self.theta_possible

        Parameters
        ----------
        q : list
            The array of Athena data to be integrated.

        Returns
        -------
        list
            The array of integrated data.
        """
        
        self.spherical_grid(get_r=True, get_theta=True, get_dphi=True)
        intq_r_sintheta_dphi = np.zeros((len(self.possible_theta), len(self.possible_r))) #these are presized without phi
        q_r_sintheta_dphi = q * self.r * np.sin(self.theta) * self.dphi
        for n in range(self.NumMeshBlocks):
            for t in range(self.theta_len):
                theta_loc = int(np.argwhere(self.possible_theta == self.theta_primitive[n, t]))
                for r in range(self.r_len):
                    r_loc = int(np.argwhere(self.possible_r == self.r_primitive[n][r]))
                    intq_r_sintheta_dphi[theta_loc, r_loc] += np.sum(q_r_sintheta_dphi[n, :, t, r])
        return intq_r_sintheta_dphi

    def _ARCHAIC_shell_integral(self, q, axis='theta', bounds=[np.pi/4,3/4*np.pi]):
        """
        Shell integrates Athena data arrays. The first integral is azmuthal then it can integrate over theta or r

        Parameters
        ----------
        q : list
            The data array to be integrated
        axis : str, default 'theta'
            The axis you wish to perform the second integration over, accepts 'theta' or 'r'
        """
        if axis == 'theta':
            az_intq = self.azmuthal_integral(q)
            shell_intq = np.zeros(len(self.possible_r))
            lower_bound = np.argmin(abs(self.possible_theta - bounds[0]))
            upper_bound = np.argmin(abs(self.possible_theta - bounds[1]))
            for r in range(len(self.possible_r)):
                shell_intq[r] = np.sum(az_intq[lower_bound:upper_bound,r] * self.possible_r[r] * self.possible_dtheta_primitive[lower_bound:upper_bound])
        if axis == 'r': #unchecked unsure if this works
            az_intq = self.azmuthal_integral(q)
            #dr still missing
            shell_intq = np.sum(az_intq, axis = 1)
        return shell_intq, az_intq

    def integrate(self, q:list, variable:str, bounds:list=None, second_variable:str=None, second_bounds:list=None, third_bounds:list=None, intermediates:bool=False, RUSTY:bool=True)->list:
        """
        Integrates athena data and outputs an array of integrated data.
        KNOWN BUG WITH CYLINDRICAL PHI INTEGRAL IN Z-PHI-R PATHWAY (fixed?)

        Parameters
        ----------
        q : list
            The data to be integrated. Accepts either Athena data array, or a single integer or float.
        variable : str
            The variable to be integrated over, accepts conventional variable names such as 'r', 'theta', 'z' or 'phi'. Also accepts the keyword "all" to perform a integral over the entire volume or "shell" to integrate through all except r.
        bounds : list of float, default None
            Sets the bounds of the integral. If bounds is None or "Full" it will integrate over the entire domain of the simulation. Bounds should be provided as an array of the form [lower-bound, upper_bound].
        second_variable : str, default None
            The second variable to be integrated over, if this field is not None a double integral will be performed over first variable and then second_variable.
        second_bounds : list of float, default None
            Sets the bounds of the second integral in the double integral, if secound_bounds is None or "Full" it will integrate over the entire domain of the simulation. Bounds should be provided as an array of the form [lower-bound, upper_bound].
        third_bounds : list of float, default None
            Sets the bounds of the third integral, if this is not None a triple integral will be performed over all 3 coordinates. Bounds should be provided as an array of the form [lower-bound, upper_bound].
            Unlike the second_bounds kwarg if third_bounds is left as None no triple integral is computed. If you would like to integrate over the entire domain use third_bounds="Full".
        intermediates : bool, default False
            Determines if intermediate integrals will be returned.

        Returns
        -------
        list
            The arrays of integrated data. If you only integrated once one array will be returned. For multiple integrals only the final result is returned, unless intermediates is True in which case the following applies:
            If you did a double integral the final result and the intermediate result from the first integral will be returned in that order,
            if you did a triple integral the final result, the intermediate result of the double integral and the intermediate result of the first integral will all be returned in that order.
        """
        self._axes()
        
        if variable == "All" or variable == "all":
            if self.gridtype == "Spherical":
                variable = "phi"
                second_variable = "theta"
                third_bounds = "Full"
            if self.gridtype == "Cylindrical":
                variable = "phi"
                second_variable = "z"
                third_bounds = "Full"
        if variable == "Shell" or variable == "shell":
            if self.gridtype == "Spherical":
                variable = "phi"
                second_variable = "theta"
                second_bounds=[np.pi/4, 3*np.pi/4]
                warning.warn("enforced truncated shell bound pi/4 to 3pi/4")
            if self.gridtype == "Cylindrical":
                variable = "phi"
                second_variable = "z"
        if variable == "Vert" or variable == "vert":
            if self.gridtype == "Spherical":
                variable = "theta"
            if self.gridtype == "Cylindrical":
                variable = "z"

        #check for improper use pt 1: bad parameters
        if variable == second_variable:
            raise Exception("You can't integrate over the same variable twice")
        if second_variable is None and second_bounds is not None:
            warnings.warn("It appears you gave bounds for a second integral but didn't tell me what variable the second integral should be integrated over. As a result Athena Analysis will not compute a second integral")
        
        if self.gridtype == "Cylindrical":
            if (variable == "r" and second_variable == "phi") or (((variable == "r" and second_variable == "z") or (variable == "z" and second_variable == "r")) and third_bounds is not None):
                raise Exception("You cannot integrate through r before phi")
        if self.gridtype == "Spherical":
            if (variable == "r" and second_variable is not None) or (second_variable == "r" and third_bounds is not None):
                raise Exception("You cannot integrate through r before theta or phi")
            if (variable == "theta" and second_variable == "phi") or ((variable == "theta" and second_variable == "r") and third_bounds is not None):
                raise Exception("You cannot integrate through theta before phi")

        if intermediates and RUSTY:
            warnings.warn("Rusty infrastructure can't handle intermediates, defaulting to python")
            RUSTY = False

        if RUSTY:
            if isinstance(q, int) or isinstance(q, float):
                q = np.full(self.array_size, q, dtype="float64")
            else:
                q = np.array(q, dtype="float64")
            
            if second_variable is not None:
                if third_bounds is not None:
                    if self.gridtype == "Cylindrical":
                        if variable == "z" or second_variable == "z":
                            variables = [variable, second_variable, "r"]
                        elif second_variable == "r":
                            variables = [variable, second_variable, "z"]
                    if self.gridtype == "Spherical":
                        variables = ["theta", "phi", "r"] # uses math convention so phi and theta are flipped
                    
                    if third_bounds == "Full":
                        third_bounds = [self.possible_vars[variables[2]][0], self.possible_vars[variables[2]][-1]]
                    if second_bounds == "Full" or second_bounds is None:
                        second_bounds = [self.possible_vars[variables[1]][0], self.possible_vars[variables[1]][-1]]
                    if bounds == "Full" or bounds is None:
                        bounds = [self.possible_vars[variables[0]][0], self.possible_vars[variables[0]][-1]]
                    bounds_tpl = ((bounds[0]-1, bounds[1]+1), (second_bounds[0]-1, second_bounds[1]+1), (third_bounds[0]-1, third_bounds[1]+1))
                else:
                    variables = [variable, second_variable]

                    if second_bounds == "Full" or second_bounds is None:
                        second_bounds = [self.possible_vars[variables[1]][0], self.possible_vars[variables[1]][-1]]
                    if bounds == "Full" or bounds is None:
                        bounds = [self.possible_vars[variables[0]][0], self.possible_vars[variables[0]][-1]]
                    bounds_tpl = ((bounds[0]-1, bounds[1]+1), (second_bounds[0]-1, second_bounds[1]+1))
            else:
                variables = [variable]

                if bounds == "Full" or bounds is None:
                    bounds = [self.possible_vars[variables[0]][0], self.possible_vars[variables[0]][-1]]
                bounds_tpl = (bounds[0]-1, bounds[1]+1)
            
            if len(variables) == 3:
                if self.gridtype == "Spherical":
                    self.native_grid(get_r=True, get_dr=True, get_phi=True, get_dphi=True, get_theta=True, get_dtheta=True)
                    return rusty_athena.volume_integral(q, self.phi_primitive, self.theta_primitive, self.r_primitive,
                                                self.dphi, self.dtheta, self.dr,
                                                self.possible_phi, self.possible_theta, self.possible_r,
                                                self.possible_dphi_primitive, self.possible_dtheta_primitive, self.possible_dr_primitive,
                                                variables, bounds_tpl)
                if self.gridtype == "Cylindrical":
                    self.native_grid(get_r=True, get_dr=True, get_phi=True, get_dphi=True, get_z=True, get_dz=True)
                    return rusty_athena.volume_integral(q, self.z_primitive, self.phi_primitive, self.r_primitive,
                                                self.dz, self.dphi, self.dr,
                                                self.possible_z, self.possible_phi, self.possible_r,
                                                self.possible_dz_primitive, self.possible_dphi_primitive, self.possible_dr_primitive,
                                                variables, bounds_tpl)
            elif len(variables) == 2:
                if self.gridtype == "Spherical":
                    self.native_grid(get_r=True, get_dr=True, get_phi=True, get_dphi=True, get_theta=True, get_dtheta=True)
                    return rusty_athena.area_integral(q, self.phi_primitive, self.theta_primitive, self.r_primitive,
                                                self.dphi, self.dtheta, self.dr,
                                                self.possible_phi, self.possible_theta, self.possible_r,
                                                self.possible_dphi_primitive, self.possible_dtheta_primitive, self.possible_dr_primitive,
                                                variables, bounds_tpl)
                if self.gridtype == "Cylindrical":
                    self.native_grid(get_r=True, get_dr=True, get_phi=True, get_dphi=True, get_z=True, get_dz=True)
                    return rusty_athena.area_integral(q, self.z_primitive, self.phi_primitive, self.r_primitive,
                                                self.dz, self.dphi, self.dr,
                                                self.possible_z, self.possible_phi, self.possible_r,
                                                self.possible_dz_primitive, self.possible_dphi_primitive, self.possible_dr_primitive,
                                                variables, bounds_tpl)
            elif len(variables) == 1:
                if self.gridtype == "Spherical":
                    self.native_grid(get_r=True, get_dr=True, get_phi=True, get_dphi=True, get_theta=True, get_dtheta=True)
                    return rusty_athena.single_integral(q, self.phi_primitive, self.theta_primitive, self.r_primitive,
                                                self.dphi, self.dtheta, self.dr,
                                                self.possible_phi, self.possible_theta, self.possible_r,
                                                variables, bounds_tpl)
                if self.gridtype == "Cylindrical":
                    self.native_grid(get_r=True, get_dr=True, get_phi=True, get_dphi=True, get_z=True, get_dz=True)
                    return rusty_athena.single_integral(q, self.z_primitive, self.phi_primitive, self.r_primitive,
                                                self.dz, self.dphi, self.dr,
                                                self.possible_z, self.possible_phi, self.possible_r,
                                                variables, bounds_tpl)
        else:
            if self.gridtype == "Spherical":
                #check if no integral has previously occured
                if isinstance(q, int) or isinstance(q, float) or q.shape == self.array_size:
                    if variable == 'r':
                        self.spherical_grid(get_theta=True, get_phi=True, get_dr=True)
                        intq_dr = np.zeros((len(self.possible_phi), len(self.possible_theta)))
                        q_dr = q * self.dr
                        if bounds != None and bounds != "Full":
                            self.get_slice_idx('r', (bounds[0]+bounds[1])/2, (bounds[1]-bounds[0])/2)
                            q_dr = q_dr[self.slice_idx]
                        for n in range(self.NumMeshBlocks):
                            for p in range(self.phi_len):
                                phi_loc = int(np.searchsorted(self.possible_phi, self.phi_primitive[n, p]))
                                for t in range(self.theta_len):
                                    theta_loc = int(np.searchsorted(self.possible_theta, self.theta_primitive[n][t]))
                                    intq_dr[phi_loc, theta_loc] += np.sum(q_dr[n, p, t, :])
                        if second_variable is None:
                            return intq_dr
                        if second_variable is not None:
                            raise Exception("You should always integrate through r last")
                    if variable == 'theta':
                        self.spherical_grid(get_r=True, get_phi=True, get_dtheta=True)
                        intq_r_dtheta = np.zeros((len(self.possible_phi), len(self.possible_r)))
                        q_r_dtheta = q * self.r * self.dtheta
                        if bounds != None and bounds != "Full":
                            self.get_slice_idx('theta', (bounds[0]+bounds[1])/2, (bounds[1]-bounds[0])/2)
                            q_r_dtheta = q_r_dtheta[self.slice_idx]
                        for n in range(self.NumMeshBlocks):
                            for p in range(self.phi_len):
                                phi_loc = int(np.searchsorted(self.possible_phi, self.phi_primitive[n, p]))
                                for r in range(self.r_len):
                                    r_loc = int(np.searchsorted(self.possible_r, self.r_primitive[n][r]))
                                    intq_r_dtheta[phi_loc, r_loc] += np.sum(q_r_dtheta[n, p, :, r])
                        if second_variable is None:
                            return intq_r_dtheta
                        if second_variable == 'phi':
                            raise Exception("For shell integrals please integrate through phi first")
                        if second_variable == 'r':
                            intq2_dr = np.zeros(len(self.possible_phi))
                            if second_bounds is None or second_bounds == "Full":
                                [lower_bound, upper_bound] = [0, -1]
                            if second_bounds is not None:
                                lower_bound = np.argmin(abs(self.possible_r - second_bounds[0]))
                                upper_bound = np.argmin(abs(self.possible_r - second_bounds[1]))
                            for p in range(len(self.possible_phi)):
                                intq2_dr[p] = np.sum(intq_r_dtheta[p,lower_bound:upper_bound] * self.possible_dr_primitive[lower_bound:upper_bound])
                            if third_bounds is None:
                                if intermediates == True:
                                    return intq2_dr, intq_r_dtheta
                                else:
                                    return intq2_dr
                            if third_bounds is not None:
                                raise Exception("r should always been the last variable integrated through, if you wish to take a triple integral please do it in phi theta r order.")
                    if variable == 'phi':
                        self.spherical_grid(get_r=True, get_theta=True, get_dphi=True)
                        intq_r_sintheta_dphi = np.zeros((len(self.possible_theta), len(self.possible_r)))
                        q_r_sintheta_dphi = q * self.r * np.sin(self.theta) * self.dphi
                        if bounds != None and bounds != "Full":
                            self.get_slice_idx('phi', (bounds[0]+bounds[1])/2, (bounds[1]-bounds[0])/2)
                            q_r_sintheta_dphi = q_r_sintheta_dphi[self.slice_idx]
                        for n in range(self.NumMeshBlocks):
                            for t in range(self.theta_len):
                                theta_loc = int(np.searchsorted(self.possible_theta, self.theta_primitive[n, t]))
                                for r in range(self.r_len):
                                    r_loc = int(np.searchsorted(self.possible_r, self.r_primitive[n][r]))
                                    intq_r_sintheta_dphi[theta_loc, r_loc] += np.sum(q_r_sintheta_dphi[n, :, t, r])
                        if second_variable is None:
                            return intq_r_sintheta_dphi
                        if second_variable == 'theta':
                            intq2_r_dtheta = np.zeros(len(self.possible_r))
                            if second_bounds is None or second_bounds == "Full":
                                [lower_bound, upper_bound] = [0, -1]
                            if second_bounds is not None:
                                lower_bound = np.argmin(abs(self.possible_theta - second_bounds[0]))
                                upper_bound = np.argmin(abs(self.possible_theta - second_bounds[1]))
                            for r in range(len(self.possible_r)):
                                intq2_r_dtheta[r] = np.sum(intq_r_sintheta_dphi[lower_bound:upper_bound,r] * self.possible_r[r] * self.possible_dtheta_primitive[lower_bound:upper_bound])
                            if third_bounds is None:
                                if intermediates == True:
                                    return intq2_r_dtheta, intq_r_sintheta_dphi
                                else:
                                    return intq2_r_dtheta
                            if third_bounds is not None:
                                intq3_dr = np.sum(intq2_r_dtheta * self.possible_dr_primitive)
                                if intermediates == True:
                                    return intq3_dr, intq2_r_dtheta, intq_r_sintheta_dphi
                                else:
                                    return intq3_dr
                        if second_variable == 'r':
                            intq2_dr = np.zeros(len(self.possible_theta))
                            if second_bounds is None or second_bounds == "Full":
                                [lower_bound, upper_bound] = [0, -1]
                            if second_bounds is not None:
                                lower_bound = np.argmin(abs(self.possible_r - second_bounds[0]))
                                upper_bound = np.argmin(abs(self.possible_r - second_bounds[1]))
                            for t in range(len(self.possible_theta)):
                                intq2_dr[t] = np.sum(intq_r_sintheta_dphi[t,lower_bound:upper_bound] * self.possible_dr_primitive[lower_bound:upper_bound])
                            if third_bounds is None:
                                if intermediates == True:
                                    return intq2_dr, intq_r_sintheta_dphi
                                else:
                                    return intq2_dr
                            if third_bounds is not None:
                                raise Exception("r should always been the last variable integrated through, if you wish to take a triple integral please do it in phi theta r order.")
                #check if q has already been integrated through phi
                if q.shape == (len(self.possible_theta), len(self.possible_r)):
                    if len(self.possible_theta) == len(self.possible_phi):
                        warnings.warn("Athena Analysis can't tell if this was previously integrated through theta or phi so it's assuming phi. If this is wrong redo both this and the original integral using the second_variable kwarg of this method to take the double integral correctly in one step")
                    if variable == 'theta':
                        intq2_r_dtheta = np.zeros(len(self.possible_r))
                        if bounds is None or bounds == "Full":
                                [lower_bound, upper_bound] = [0, -1]
                        if bounds is not None:
                            lower_bound = np.argmin(abs(self.possible_theta - bounds[0]))
                            upper_bound = np.argmin(abs(self.possible_theta - bounds[1]))
                        for r in range(len(self.possible_r)):
                            intq2_r_dtheta[r] = np.sum(q[lower_bound:upper_bound,r] * self.possible_r[r] * self.possible_dtheta_primitive[lower_bound:upper_bound])
                        if second_variable is None:
                            return intq2_r_dtheta
                        if second_varibale == 'r':
                            if second_bounds == "Full" or second_bounds is None:
                                [lower_bound, upper_bound] = [0, -1]
                            else:
                                lower_bound = np.argmin(abs(self.possible_r - second_bounds[0]))
                                upper_bound = np.argmin(abs(self.possible_r - second_bounds[1]))
                            intq3_dr = np.sum(q * self.possible_dr_primitive)
                            if third_bounds is None:
                                if intermediates == True:
                                    return intq3_dr, intq2_r_dtheta
                                else:
                                    return intq3_dr
                            if third_bounds is not None:
                                raise Exception("r should always been the last variable integrated through, if you wish to take a triple integral please do it in phi theta r order.")
                        if second_variable == 'phi':
                            raise Exception("phi should be integrated through before theta, if you wish to take a double or triple integral please do it in phi theta r order.")
                    if variable == 'r':
                        intq2_dr = np.zeros(len(self.possible_theta))
                        if bounds is None or bounds == "Full":
                                [lower_bound, upper_bound] = [0, -1]
                        if bounds is not None:
                            lower_bound = np.argmin(abs(self.possible_r - bounds[0]))
                            upper_bound = np.argmin(abs(self.possible_r - bounds[1]))
                        for t in range(len(self.possible_theta)):
                            intq2_dr[t] = np.sum(q[t,lower_bound:upper_bound] * self.possible_dr_primitive[lower_bound:upper_bound])
                        if second_variable is None:
                            return intq2_dr
                        if second_variabel is not None:
                            raise Exception("r should always been the last variable integrated through, if you wish to take a triple integral please do it in phi theta r order.")
                    if variable == 'phi':
                        raise Exception("I think you've already integrated through this variable, you can't integrate over the same axis twice. If I'm wrong use the inbuilt multiple integral functionality to ensure your integral is taken correctly")
                #check if q has already been integrated through theta
                if q.shape == (len(self.possible_phi), len(self.possible_r)):
                    if len(self.possible_theta) == len(self.possible_phi):
                        warnings.warn("Athena Analysis can't tell if this was previously integrated through theta or phi so it's assuming theta. If this is wrong redo both this and the original integral using the second_variable kwarg of this method to take the double integral correctly in one step")
                    if len(self.possible_theta) == len(self.possible_r):
                        warnings.warn("Athena Analysis can't tell if this was previously integrated through theta or r so it's assuming theta. If this is wrong redo both this and the original integral using the second_variable kwarg of this method to take the double integral correctly in one step")
                    if variable == 'phi':
                        raise Exception("For shell integrals please integrate through phi first")
                    if variable == 'r':
                        intq2_dr = np.zeros(len(self.possible_phi))
                        if bounds is None or bounds == "Full":
                                [lower_bound, upper_bound] = [0, -1]
                        if bounds is not None:
                            lower_bound = np.argmin(abs(self.possible_r - bounds[0]))
                            upper_bound = np.argmin(abs(self.possible_r - bounds[1]))
                        for p in range(len(self.possible_phi)):
                            intq2_dr[p] = np.sum(q[t,lower_bound:upper_bound] * self.possible_dr_primitive[lower_bound:upper_bound])
                        return intq2_dr
                    if variable == 'phi':
                        raise Exception("phi should be integrated through before theta, if you wish to take a double or triple integral please do it in phi theta r order.")
                    if variable == 'theta':
                        raise Exception("I think you've already integrated through this variable, you can't integrate over the same axis twice. If I'm wrong use the inbuilt multiple integral functionality to ensure your integral is taken correctly")
                #check if q has already been integrated through r
                if q.shape == (len(self.possible_phi), len(self.possible_thets)):
                    raise Exception("Anthena Analysis thinks this data has already been integrated through r and r should always been the last variable integrated through. If this is wrong redo both this and the original integral using the second_variable kwarg of this method to take the double integral correctly in one step")
                #check if q has already been integrated through phi and theta
                if q.shape == (len(self.possible_r)):
                    if len(self.possible_r) == len(self.possible_theta):
                        warnings.warn("Athena Analysis can't tell if this was previously integrated through phi and theta or phi and r so it's assuming r. If this is wrong redo both this and the original integrals using the inbuilt triple integral capability of this method to take the triple integral correctly in one step")
                    if len(self.possible_r) == len(self.possible_phi):
                        warnings.warn("Athena Analysis can't tell if this was previously integrated through phi and theta or theta and r so it's assuming r. If this is wrong redo both this and the original integrals using the inbuilt triple integral capability of this method to take the triple integral correctly in one step")
                    if variable == 'r':
                        if bounds == "Full" or bounds is None:
                                    [lower_bound, upper_bound] = [0, -1]
                        else:
                            lower_bound = np.argmin(abs(self.possible_r - third_bounds[0]))
                            upper_bound = np.argmin(abs(self.possible_r - third_bounds[1]))
                        intq3_dr = np.sum(q * self.possible_dr_primitive)
                        return intq3_dr
                    if variable != 'r':
                        raise Exception("I think you've already integrated through this variable, you can't integrate over the same axis twice. If I'm wrong use the inbuilt multiple integral functionality to ensure your integral is taken correctly")
                #check for improper use pt 2: bad triple integrals
                if q.shape == (len(self.possible_theta)):
                    raise Exception("It appears you trying to finish a triple integral with theta as the final variable integrated over. r must be the final variable integrated over. It is recommended you redo both this and the original integrals using the inbuilt triple integral capability of this method to take the triple integral correctly in one step")
                if q.shape == (len(self.possible_phi)):
                    raise Exception("It appears you trying to finish a triple integral with phi as the final variable integrated over. r must be the final variable integrated over. It is recommended you redo both this and the original integrals using the inbuilt triple integral capability of this method to take the triple integral correctly in one step")                
            if self.gridtype == "Cylindrical":
                if isinstance(q, int) or isinstance(q, float) or q.shape == self.array_size:
                    if variable == 'phi':
                        self.native_grid(get_r=True, get_dphi=True)
                        intq_r_dphi = np.zeros((len(self.possible_z), len(self.possible_r)))
                        q_r_dphi = q * self.r * self.dphi
                        if bounds != None and bounds != "Full":
                            self.get_slice_idx('phi', (bounds[0]+bounds[1])/2, (bounds[1]-bounds[0])/2)
                            q_r_dphi = q_r_dphi[self.slice_idx]
                        for n in range(self.NumMeshBlocks):
                            for z in range(self.z_len):
                                z_loc = int(np.searchsorted(self.possible_z, self.z_primitive[n, z]))
                                for r in range(self.r_len):
                                    r_loc = int(np.searchsorted(self.possible_r, self.r_primitive[n][r]))
                                    intq_r_dphi[z_loc, r_loc] += np.sum(q_r_dphi[n, z, :, r])
                        if second_variable is None:
                            return intq_r_dphi
                        if second_variable == 'z':
                            intq2_dz = np.zeros(len(self.possible_r))
                            if second_bounds is None or second_bounds == "Full":
                                [lower_bound, upper_bound] = [0, -1]
                            if second_bounds is not None:
                                lower_bound = np.argmin(abs(self.possible_z - second_bounds[0]))
                                upper_bound = np.argmin(abs(self.possible_z - second_bounds[1]))
                            for r in range(len(self.possible_r)):
                                intq2_dz[r] = np.sum(intq_r_dphi[lower_bound:upper_bound,r] * self.possible_dz_primitive[lower_bound:upper_bound])
                            if third_bounds is None:
                                if intermediates == True:
                                    return intq2_dz, intq_r_dphi
                                else:
                                    return intq2_dz
                            if third_bounds is not None:
                                if third_bounds == "Full":
                                    [lower_bound, upper_bound] = [0, -1]
                                if third_bounds != "Full":
                                    lower_bound = np.argmin(abs(self.possible_r - third_bounds[0]))
                                    upper_bound = np.argmin(abs(self.possible_r - third_bounds[1]))
                                intq3_dr = np.sum(intq2_dz[lower_bound:upper_bound] * self.possible_dr_primitive[lower_bound:upper_bound])
                                if intermediates == True:
                                    return intq3_dr, intq2_dz, intq_r_dphi
                                else:
                                    return intq3_dr
                        if second_variable == 'r':
                            intq2_dr = np.zeros(len(self.possible_z))
                            if second_bounds is None or second_bounds == "Full":
                                [lower_bound, upper_bound] = [0, -1]
                            if second_bounds is not None:
                                lower_bound = np.argmin(abs(self.possible_r - second_bounds[0]))
                                upper_bound = np.argmin(abs(self.possible_r - second_bounds[1]))
                            for z in range(len(self.possible_z)):
                                intq2_dr[z] = np.sum(intq_r_dphi[z,lower_bound:upper_bound] * self.possible_dr_primitive[lower_bound:upper_bound])
                            if third_bounds is None:
                                if intermediates == True:
                                    return intq2_dr, intq_r_dphi
                                else:
                                    return intq2_dr
                            if third_bounds is not None:
                                if third_bounds == "Full":
                                    [lower_bound, upper_bound] = [0, -1]
                                if third_bounds != "Full":
                                    lower_bound = np.argmin(abs(self.possible_z - second_bounds[0]))
                                    upper_bound = np.argmin(abs(self.possible_z - second_bounds[1]))
                                intq3_dz = np.sum(intq2_dr[lower_bound:upper_bound] * self.possible_dz_primitive[lower_bound:upper_bound])
                                if intermediates == True:
                                    return intq3_dz, intq2_dr, intq_r_dphi
                                else:
                                    return intq3_dz
                    if variable == 'z':
                        self.native_grid(get_dz=True)
                        intq_dz = np.zeros((len(self.possible_phi), len(self.possible_r)))
                        q_dz = q * self.dz
                        if bounds != None and bounds != "Full":
                            self.get_slice_idx('z', (bounds[0]+bounds[1])/2, (bounds[1]-bounds[0])/2)
                            q_dz = q_dz[self.slice_idx]
                        for n in range(self.NumMeshBlocks):
                            for p in range(self.phi_len):
                                phi_loc = int(np.searchsorted(self.possible_phi, self.phi_primitive[n, p]))
                                for r in range(self.r_len):
                                    r_loc = int(np.searchsorted(self.possible_r, self.r_primitive[n][r]))
                                    intq_dz[phi_loc, r_loc] += np.sum(q_dz[n, :, p, r])
                        if second_variable is None:
                            return intq_dz
                        if second_variable == 'phi':
                            intq2_r_dphi = np.zeros(len(self.possible_r))
                            if second_bounds is None or second_bounds == "Full":
                                [lower_bound, upper_bound] = [0, -1]
                            if second_bounds is not None:
                                lower_bound = np.argmin(abs(self.possible_phi - second_bounds[0]))
                                upper_bound = np.argmin(abs(self.possible_phi - second_bounds[1]))
                            for r in range(len(self.possible_r)):
                                intq2_r_dphi[r] = np.sum(intq_dz[lower_bound:upper_bound,r] * self.possible_r[r] * self.possible_dphi_primitive[lower_bound:upper_bound])
                            if third_bounds is None:
                                if intermediates == True:
                                    return intq2_r_dphi, intq_dz
                                else:
                                    return intq2_r_dphi
                            if third_bounds is not None:
                                if third_bounds == "Full":
                                    [lower_bound, upper_bound] = [0, -1]
                                if third_bounds != "Full":
                                    lower_bound = np.argmin(abs(self.possible_r - second_bounds[0]))
                                    upper_bound = np.argmin(abs(self.possible_r - second_bounds[1]))
                                intq3_dr = np.sum(intq2_r_dphi[lower_bound:upper_bound] * self.possible_dr_primitive[lower_bound:upper_bound])
                                if intermediates == True:
                                    return intq3_dr, intq2_r_dphi, intq_dz
                                else:
                                    return intq3_dr
                        if second_variable == 'r':
                            intq2_dr = np.zeros(len(self.possible_phi))
                            if second_bounds is None or second_bounds == "Full":
                                [lower_bound, upper_bound] = [0, -1]
                            if second_bounds is not None:
                                lower_bound = np.argmin(abs(self.possible_r - second_bounds[0]))
                                upper_bound = np.argmin(abs(self.possible_r - second_bounds[1]))
                            for p in range(len(self.possible_p)):
                                intq2_dr[p] = np.sum(intq_dz[p,lower_bound:upper_bound] * self.possible_dr_primitive[lower_bound:upper_bound])
                            if third_bounds is None:
                                if intermediates == True:
                                    return intq2_dr, intq_dz
                                else:
                                    return intq2_dr
                            if third_bounds is not None:
                                raise Exception("r must integrated through after phi, please reorder your integral (you are currently attempting dphidrdz order).")
                    if variable == 'r':
                        self.native_grid(get_dr=True)
                        intq_dr = np.zeros((len(self.possible_z), len(self.possible_phi)))
                        q_dr = q * self.dr
                        if bounds != None and bounds != "Full":
                            self.get_slice_idx('r', (bounds[0]+bounds[1])/2, (bounds[1]-bounds[0])/2)
                            q_dr = q_dr[self.slice_idx]
                        for n in range(self.NumMeshBlocks):
                            for z in range(self.z_len):
                                z_loc = int(np.searchsorted(self.possible_z, self.z_primitive[n, z]))
                                for p in range(self.phi_len):
                                    phi_loc = int(np.searchsorted(self.possible_phi, self.phi_primitive[n][p]))
                                    intq_dr[z_loc, phi_loc] += np.sum(q_dr[n, z, p, :])
                        if second_variable is None:
                            return intq_dr
                        if second_variable == 'z':
                            intq2_dz = np.zeros(len(self.possible_phi))
                            if second_bounds is None or second_bounds == "Full":
                                [lower_bound, upper_bound] = [0, -1]
                            if second_bounds is not None:
                                lower_bound = np.argmin(abs(self.possible_z - second_bounds[0]))
                                upper_bound = np.argmin(abs(self.possible_z - second_bounds[1]))
                            for p in range(len(self.possible_phi)):
                                intq2_dz[p] = np.sum(intq_dr[lower_bound:upper_bound,p] * self.possible_dz_primitive[lower_bound:upper_bound])
                            if third_bounds is None:
                                if intermediates == True:
                                    return intq2_dz, intq_dr
                                else:
                                    return intq2_dz
                            if third_bounds is not None:
                                raise Exception("r must integrated through after phi, please reorder your integral (you are currently attempting drdzdphi order).")
                        if second_variable == 'phi':
                            raise Exception("r must integrated through after phi, please reorder your integral (you are currently attempting drdphi order).")
                #check if q has already been integrated through z
                if q.shape == (len(self.possible_phi), len(self.possible_r)):
                    if len(self.possible_phi) == len(self.possible_z):
                        warnings.warn("Athena Analysis can't tell if this was previously integrated through z or phi so it's assuming z. If this is wrong redo both this and the original integral using the second_variable kwarg of this method to take the double integral correctly in one step")
                    if variable == 'phi':
                        intq2_r_dphi = np.zeros(len(self.possible_r))
                        if bounds is None or bounds == "Full":
                                [lower_bound, upper_bound] = [0, -1]
                        if bounds is not None:
                            lower_bound = np.argmin(abs(self.possible_phi - bounds[0]))
                            upper_bound = np.argmin(abs(self.possible_phi - bounds[1]))
                        for r in range(len(self.possible_r)):
                            intq2_r_dphi[r] = np.sum(q[lower_bound:upper_bound,r] * self.possible_r[r] * self.possible_dphi_primitive[lower_bound:upper_bound])
                        if second_variable is None:
                            return intq2_r_dphi
                        if second_variable == 'r':
                            if second_bounds == "Full" or second_bounds is None:
                                [lower_bound, upper_bound] = [0, -1]
                            else:
                                lower_bound = np.argmin(abs(self.possible_r - second_bounds[0]))
                                upper_bound = np.argmin(abs(self.possible_r - second_bounds[1]))
                            intq3_dr = np.sum(q * self.possible_dr_primitive)
                            if third_bounds is None:
                                if intermediates == True:
                                    return intq3_dr, intq2_r_dphi
                                else:
                                    return intq3_dr
                            if third_bounds is not None:
                                raise Exception("hmm... I think you've already integrated through every axis, not sure what you're trying to do here")
                    if variable == 'r':
                        intq2_dr = np.zeros(len(self.possible_phi))
                        if bounds is None or bounds == "Full":
                                [lower_bound, upper_bound] = [0, -1]
                        if bounds is not None:
                            lower_bound = np.argmin(abs(self.possible_r - bounds[0]))
                            upper_bound = np.argmin(abs(self.possible_r - bounds[1]))
                        for p in range(len(self.possible_phi)):
                            intq2_dr[p] = np.sum(q[p,lower_bound:upper_bound] * self.possible_dr_primitive[lower_bound:upper_bound])
                        if second_variable is None:
                            return intq2_dr
                        if second_variable == 'phi':
                            raise Exception("r must integrated through after phi, please reorder your integral (you are currently attempting dzdrdphi order).")
                        if second_variable == 'z':
                            raise Exception("I think you've already integrated through this variable, you can't integrate over the same axis twice. If I'm wrong use the inbuilt multiple integral functionality to ensure your integral is taken correctly")
                    if variable == 'z':
                        raise Exception("I think you've already integrated through this variable, you can't integrate over the same axis twice. If I'm wrong use the inbuilt multiple integral functionality to ensure your integral is taken correctly")
                #check if q has already been integrated through phi
                if q.shape == (len(self.possible_z), len(self.possible_r)):
                    if len(self.possible_r) == len(self.possible_phi):
                        warnings.warn("Athena Analysis can't tell if this was previously integrated through phi or r so it's assuming phi. If this is wrong redo both this and the original integral using the second_variable kwarg of this method to take the double integral correctly in one step")
                    if variable == 'z':
                        intq2_dz = np.zeros(len(self.possible_z))
                        if bounds is None or bounds == "Full":
                                [lower_bound, upper_bound] = [0, -1]
                        if bounds is not None:
                            lower_bound = np.argmin(abs(self.possible_z - bounds[0]))
                            upper_bound = np.argmin(abs(self.possible_z - bounds[1]))
                        for r in range(len(self.possible_r)):
                            intq2_dz[r] = np.sum(q[lower_bound:upper_bound,r] * self.possible_dz_primitive[lower_bound:upper_bound])
                        if second_variable is None:
                            return intq2_dz
                        if second_variable == 'r':
                            if second_bounds == "Full" or second_bounds is None:
                                [lower_bound, upper_bound] = [0, -1]
                            else:
                                lower_bound = np.argmin(abs(self.possible_r - second_bounds[0]))
                                upper_bound = np.argmin(abs(self.possible_r - second_bounds[1]))
                            intq3_dr = np.sum(q * self.possible_dr_primitive)
                            if third_bounds is None:
                                if intermediates == True:
                                    return intq3_dr, intq2_dz
                                else:
                                    return intq3_dr
                            if third_bounds is not None:
                                raise Exception("hmm... I think you've already integrated through every axis, not sure what you're trying to do here")
                    if variable == 'r':
                        intq2_dr = np.zeros(len(self.possible_z))
                        if bounds is None or bounds == "Full":
                                [lower_bound, upper_bound] = [0, -1]
                        if bounds is not None:
                            lower_bound = np.argmin(abs(self.possible_r - bounds[0]))
                            upper_bound = np.argmin(abs(self.possible_r - bounds[1]))
                        for z in range(len(self.possible_z)):
                            intq2_dr[z] = np.sum(q[z,lower_bound:upper_bound] * self.possible_dr_primitive[lower_bound:upper_bound])
                        if second_variable is None:
                            return intq2_dr
                        if second_variable == 'z':
                            if second_bounds == "Full" or second_bounds is None:
                                [lower_bound, upper_bound] = [0, -1]
                            else:
                                lower_bound = np.argmin(abs(self.possible_z - second_bounds[0]))
                                upper_bound = np.argmin(abs(self.possible_z - second_bounds[1]))
                            intq3_dz = np.sum(q * self.possible_dz_primitive)
                            if third_bounds is None:
                                if intermediates == True:
                                    return intq3_dz, intq2_dr
                                else:
                                    return intq3_dr
                            if third_bounds is not None:
                                raise Exception("hmm... I think you've already integrated through every axis, not sure what you're trying to do here")
                    if variable == 'phi':
                        raise Exception("I think you've already integrated through this variable, you can't integrate over the same axis twice. If I'm wrong use the inbuilt multiple integral functionality to ensure your integral is taken correctly")
                #check if q has already been integrated through r
                if q.shape == (len(self.possible_z), len(self.possible_phi)):
                    if len(self.possible_phi) == len(self.possible_r):
                        warnings.warn("Athena Analysis can't tell if this was previously integrated through r or phi so it's assuming r. If this is wrong redo both this and the original integral using the second_variable kwarg of this method to take the double integral correctly in one step")
                    if variable == 'phi':
                        raise Exception("r must integrated through after phi, please reorder your integral (you are currently attempting drdphi order).")
                    if variable == 'z':
                        intq2_dz = np.zeros(len(self.possible_phi))
                        if bounds is None or bounds == "Full":
                                [lower_bound, upper_bound] = [0, -1]
                        if bounds is not None:
                            lower_bound = np.argmin(abs(self.possible_z - bounds[0]))
                            upper_bound = np.argmin(abs(self.possible_z - bounds[1]))
                        for p in range(len(self.possible_phi)):
                            intq2_dr[p] = np.sum(q[p,lower_bound:upper_bound] * self.possible_dz_primitive[lower_bound:upper_bound])
                        if second_variable is None:
                            return intq2_dz
                        if second_variable == 'phi':
                            raise Exception("r must integrated through after phi, please reorder your integral (you are currently attempting dzdrdphi order).")
                        if second_variable == 'r':
                            raise Exception("I think you've already integrated through this variable, you can't integrate over the same axis twice. If I'm wrong use the inbuilt multiple integral functionality to ensure your integral is taken correctly")
                    if variable == 'z':
                        raise Exception("I think you've already integrated through this variable, you can't integrate over the same axis twice. If I'm wrong use the inbuilt multiple integral functionality to ensure your integral is taken correctly")
                #check if q has already been integrated through phi and z
                if q.shape == self.possible_r.shape:
                    if len(self.possible_r) == len(self.possible_z):
                        warnings.warn("Athena Analysis can't tell if this was previously integrated through phi and z or phi and r so it's assuming phi and z. If this is wrong redo both this and the original integrals using the inbuilt triple integral capability of this method to take the triple integral correctly in one step")
                    if len(self.possible_r) == len(self.possible_phi):
                        warnings.warn("Athena Analysis can't tell if this was previously integrated through phi and z or z and r so it's assuming phi and z. If this is wrong redo both this and the original integrals using the inbuilt triple integral capability of this method to take the triple integral correctly in one step")
                    if variable == 'r':
                        if bounds == "Full" or bounds is None:
                                    [lower_bound, upper_bound] = [0, -1]
                        else:
                            lower_bound = np.argmin(abs(self.possible_r - third_bounds[0]))
                            upper_bound = np.argmin(abs(self.possible_r - third_bounds[1]))
                        intq3_dr = np.sum(q * self.possible_dr_primitive)
                        return intq3_dr
                    if variable != 'r':
                        raise Exception("I think you've already integrated through this variable, you can't integrate over the same axis twice. If I'm wrong use the inbuilt multiple integral functionality to ensure your integral is taken correctly")
                #check if q has already been integrated through phi and r
                if q.shape == self.possible_z.shape:
                    if len(self.possible_r) == len(self.possible_z):
                        warnings.warn("Athena Analysis can't tell if this was previously integrated through phi and r or z and r so it's assuming phi and r. If this is wrong redo both this and the original integrals using the inbuilt triple integral capability of this method to take the triple integral correctly in one step")
                    if len(self.possible_z) == len(self.possible_phi):
                        warnings.warn("Athena Analysis can't tell if this was previously integrated through phi and r or z and r so it's assuming phi and r. If this is wrong redo both this and the original integrals using the inbuilt triple integral capability of this method to take the triple integral correctly in one step")
                    if variable == 'z':
                        if bounds == "Full" or bounds is None:
                                    [lower_bound, upper_bound] = [0, -1]
                        else:
                            lower_bound = np.argmin(abs(self.possible_z - third_bounds[0]))
                            upper_bound = np.argmin(abs(self.possible_z - third_bounds[1]))
                        intq3_dz = np.sum(q * self.possible_dz_primitive)
                        return intq3_dz
                    if variable != 'z':
                        raise Exception("I think you've already integrated through this variable, you can't integrate over the same axis twice. If I'm wrong use the inbuilt multiple integral functionality to ensure your integral is taken correctly")
                #check for improper use pt 2: bad triple integrals
                if q.shape == self.possible_phi.shape:
                    raise Exception("It appears you trying to finish a triple integral with phi as the final variable integrated over. Phi must be integrated over before r. It is recommended you redo both this and the original integrals using the inbuilt triple integral capability of this method to take the triple integral correctly in one step")
            #raise for improper use final pt: bad format
                raise Exception("Sorry Athen Analysis doesn't recognize the format of your data. Please make sure its either an Athena data array, an array I created from this integrator, an integer, or a float. Other arrays of custom shape will not work ;(")
    
    def vec_vol_integrate(self, v, radial_bounds=None):
        """
        Integrates every component of a vector over all space then return a vector of the three integrated components
        """
        [a, b, c] = v
        if radial_bounds is not None:
            return np.array([self.integrate(a, variable='phi', second_variable='z', third_bounds=radial_bounds), self.integrate(b, variable='phi', second_variable='z', third_bounds=radial_bounds), self.integrate(c, variable='phi', second_variable='z', third_bounds=radial_bounds)])

        return np.array([self.integrate(a, variable='All'), self.integrate(b, variable='All'), self.integrate(c, variable='All')])
    
    def vec_shell_integrate(self, v):
        """
        Integrates every component of a vector over all space then return a vector of the three integrated components
        """
        [a, b, c] = v
        return np.array([self.integrate(a, variable='shell'), self.integrate(b, variable='shell'), self.integrate(c, variable='shell')])
    
    def differentiate(self, q, variable):
        """
        Differentiates the data with respect to the given variable.
        WARNING: The derivative at meshblock edges is linearly approximated, this may cause problems if the dervivative varies sharply or you take many derivatives in succession.
        
        Parameters
        ----------
        q : list
            The data to be differentiated, must be an Athena data array.
        variable : str
            The variable with which to differentiate, accepts 'r', 'theta', or 'phi'
        """        
        if variable == 'phi':
            self.native_grid(get_phi=True)
            dqdphi = np.zeros(self.array_size)
            for n in range(self.NumMeshBlocks):
                if self.gridtype=='Spherical':
                    dqdphi[n,1:,:,:] = (q[n,1:,:,:] - q[n,:-1,:,:]) / (self.phi[n,1:,:,:] - self.phi[n,:-1,:,:])
                    dqdphi[n,0,:,:] = dqdphi[n,1,:,:] - ((dqdphi[n,2,:,:]-dqdphi[n,1,:,:]) * ((self.phi[n,1,:,:]-self.phi[n,0,:,:])/(self.phi[n,2,:,:]-self.phi[n,1,:,:])))
                if self.gridtype=='Cylindrical':
                    dqdphi[n,:,1:,:] = (q[n,:,1:,:] - q[n,:,:-1,:]) / (self.phi[n,:,1:,:] - self.phi[n,:,:-1,:])
                    dqdphi[n,:,0,:] = dqdphi[n,:,1,:] - ((dqdphi[n,:,2,:]-dqdphi[n,:,1,:]) * ((self.phi[n,:,1,:]-self.phi[n,:,0,:])/(self.phi[n,:,2,:]-self.phi[n,:,1,:])))
            return dqdphi
        if variable == 'theta':
            self.native_grid(get_theta=True)
            dqdtheta = np.zeros(self.array_size)
            for n in range(self.NumMeshBlocks):
                dqdtheta[n,:,1:,:] = (q[n,:,1:,:] - q[n,:,:-1,:]) / (self.theta[n,:,1:,:] - self.theta[n,:,:-1,:])
                dqdtheta[n,:,0,:] = dqdtheta[n,:,1,:] - ((dqdtheta[n,:,2,:]-dqdtheta[n,:,1,:]) * ((self.theta[n,:,1,:]-self.theta[n,:,0,:])/(self.theta[n,:,2,:]-self.theta[n,:,1,:])))
            return dqdtheta
        if variable == 'r':
            self.native_grid(get_r=True)
            dqdr = np.zeros(self.array_size)
            for n in range(self.NumMeshBlocks):
                dqdr[n,:,:,1:] = (q[n,:,:,1:] - q[n,:,:,:-1]) / (self.r[n,:,:,1:] - self.r[n,:,:,:-1])
                dqdr[n,:,:,0] = dqdr[n,:,:,1] - ((dqdr[n,:,:,2]-dqdr[n,:,:,1]) * ((self.r[n,:,:,1]-self.r[n,:,:,0])/(self.r[n,:,:,2]-self.r[n,:,:,1])))
            return dqdr
        if variable == 'z':
            self.native_grid(get_phi=True)
            dqdz = np.zeros(self.array_size)
            for n in range(self.NumMeshBlocks):
                dqdz[n,1:,:,:] = (q[n,1:,:,:] - q[n,:-1,:,:]) / (self.z[n,1:,:,:] - self.z[n,:-1,:,:])
                dqdz[n,0,:,:] = dqdz[n,1,:,:] - ((dqdz[n,2,:,:]-dqdz[n,1,:,:]) * ((self.z[n,1,:,:]-self.z[n,0,:,:])/(self.z[n,2,:,:]-self.z[n,1,:,:])))
            return dqdz

    def get_potentials(self, get_wd_grav=False, get_companion_grav=False, get_accel=False):
        """
        Gets various potential functions

        Parameters
        ----------
        get_wd_grav : bool, default False
            Determines if white dwarf gravitational potential is calculated.
        get_companion_grav : bool, default False
            Determines if companion gravitational potential is calculated.
        get_accel : bool, default False
            Determines if ficticious potential due to acceleration is calculated.
        """
        if get_wd_grav == True and self.wd_grav_pot is None:
            self.native_grid(get_r=True)
            self.wd_grav_pot = -1 * sim.gm1 / self.r
        if get_companion_grav == True and self.companion_grav_pot is None:
            self.cart_grid(projection='3d')
            if self.NO_Z_GRAVITY:
                dist = np.sqrt((self.x - sim.bin_sep)**2 + (self.y)**2)
            else:
                print("wtf1")
                dist = np.sqrt((self.x - sim.bin_sep)**2 + (self.y)**2 + (self.z)**2)
            self.companion_grav_pot = -1 * sim.gm2 / dist
        if get_accel == True and self.accel_pot is None:
            self.cart_grid(projection='2d')
            self.accel_pot = sim.gm2 * self.x / (sim.bin_sep)**2

    def gradient(self, q, coordinates=None):
        """
        Computes the gradient of a scalar

        Parameters
        ----------
        q : list
            The data to be differentiated, must be an Athena data array.
        
        Returns
        -------
        list
            A 5 dimensional array that indexs through gradient component, then meshblock and coords as normal.
        """
        if coordinates is None:
            coordinates = self.gridtype
        if coordinates == 'Spherical' or coordinates == 'spherical':
            self.spherical_grid(get_r=True, get_theta=True)
            gradq_phi = (self.differentiate(q, 'phi')) / (self.r * np.sin(self.theta))
            gradq_theta = (self.differentiate(q, 'theta')) / self.r
            gradq_r = self.differentiate(q, 'r')
            return np.array([gradq_phi, gradq_theta, gradq_r])
        if coordinates == 'Cylindrical' or coordinates == 'cylindrical':
            self.native_grid(get_r=True)
            gradq_z = self.differentiate(q, 'z')
            gradq_phi = (self.differentiate(q, 'phi')) / (self.r)
            gradq_r = self.differentiate(q, 'r')
            return np.array([gradq_z, gradq_phi, gradq_r])
        if coordinates == 'Cartesian' or coordinates == 'cartesian':
            raise Exception("use spherical/cylindrical lol cartesian is borked")

    def divergence(self, q, coordinates=None):
        """
        Computes the divergence of a vector

        Parameters
        ----------
        q : list
            The data to be differentiated, must be an Athena data array.
        
        Returns
        -------
        float
            The divergence
        """
        if coordinates is None:
            coordinates = self.gridtype
        if coordinates == 'Spherical' or coordinates == 'spherical':
            [qp, qt, qr] = q
            self.spherical_grid(get_r=True, get_theta=True)
            divq_phi = (self.differentiate(qp, 'phi')) / (self.r * np.sin(self.theta))
            divq_theta = (self.differentiate(qt*np.sin(self.theta), 'theta')) / (self.r * np.sin(self.theta))
            divq_r = self.differentiate(qr*(aa.r**2), 'r') / (aa.r**2)
            return divq_phi + divq_theta + divq_r
        if coordinates == 'Cylindrical' or coordinates == 'cylindrical':
            [qz, qp, qr] = q
            self.native_grid(get_r=True)
            divq_z = self.differentiate(qz, 'z')
            divq_phi = (self.differentiate(qp, 'phi')) / (self.r)
            divq_r = self.differentiate(qr*self.r, 'r') / (self.r)
            return divq_z + divq_phi + divq_r
        if coordinates == 'Cartesian' or coordinates == 'cartesian':
            raise Exception("use spherical/cylindrical lol cartesian is borked")

    def tensor_divergence(self, T, coordinates=None):
        """
        Computes the divergence of a tensor

        Parameters
        ----------
        q : list of list
            The data to be differentiated, must be an Athena data array.
        
        Returns
        -------
        list
            A 4 dimensional array that contains the divergence, with indices of meshblock and coords as normal.
        """
        if coordinates is None:
            coordinates = self.gridtype
        if coordinates == 'Spherical' or coordinates == 'spherical':
            raise("code not written")
            return T[0]
        if coordinates == 'Cylindrical' or coordinates == 'cylindrical':
            [[Tzz, Tzp, Tzr], 
            [Tpz, Tpp, Tpr],
            [Trz, Trp, Trr]] = T
            self.native_grid(get_r=True)
            divT_z = self.differentiate(Trz, 'r') + (self.differentiate(Tpz, 'phi') + Trz)/self.r + self.differentiate(Tzz, 'z')
            divT_phi = self.differentiate(Trp, 'r') + (self.differentiate(Tpp, 'phi') + Trp + Tpr)/self.r + self.differentiate(Tzp, 'z')
            divT_r = self.differentiate(Trr, 'r') + (self.differentiate(Tpr, 'phi') + Trr - Tpp)/self.r + self.differentiate(Tzr, 'z')
            return np.array([divT_z, divT_phi, divT_r])
        if coordinates == 'Cartesian' or coordinates == 'cartesian':
            raise Exception("use spherical/cylindrical lol cartesian is borked")

    def material_derivative(self, a, b, coordinates=None):
        """
        Computes the spacial part of a material derivative: (a・∇)・b

        Parameters
        ----------
        b : list
            The data to be differentiated, may be a 5d array consisting of an array of 3 Athena data arrays or a normal 4d Athena data array.
        a : list
            The vector to be used as the prefered direction in the directional derivative. Must be a 5d array consisting of an array of 3 Athena data arrays.
        
        Returns
        -------
        list
            A 5 dimensional array that indexs through gradient component, then meshblock and coords as normal or a normal 4d Athena data array depending on if b is a vector or scalar.
        """
        if coordinates is None:
            coordinates = self.gridtype
        b = np.array(b)
        if b.shape == self.vector_array_size:
            if coordinates == 'Spherical' or coordinates == 'spherical':
                matderiv_phi = ((a[0]*b[1]/(np.tan(self.theta)*self.r))
                    + a[0]*b[2]/self.r
                    + a[0] * (self.differentiate(b[0], 'phi') / (self.r * np.sin(self.theta)))
                    + a[1] * (self.differentiate(b[0], 'theta') / self.r)
                    + a[2] * self.differentiate(b[0], 'r'))
                matderiv_theta = ((-1*a[0]*b[0]/(np.tan(self.theta)*self.r))
                    + (a[1]*b[2]/self.r)
                    + a[0] * (self.differentiate(b[1], 'phi') / (self.r * np.sin(self.theta)))
                    + a[1] * (self.differentiate(b[1], 'theta') / self.r)
                    + a[2] * self.differentiate(b[1], 'r'))
                matderiv_r = ((-1*(a[1]*b[1]+a[0]*b[0])/self.r)
                    + a[0] * (self.differentiate(b[2], 'phi') / (self.r * np.sin(self.theta)))
                    + a[1] * (self.differentiate(b[2], 'theta') / self.r)
                    + a[2] * self.differentiate(b[2], 'r'))
                return np.array([matderiv_phi, matderiv_theta, matderiv_r])
            if coordinates == 'Cylindrical' or coordinates == 'cylindrical':
                matderiv_z = (a[0] * self.differentiate(b[0], 'z') 
                    + (a[1]/self.r) * self.differentiate(b[0], 'phi') 
                    + a[2] * self.differentiate(b[0], 'r'))
                matderiv_phi = ((a[1]*b[2]/self.r)
                    + a[0] * self.differentiate(b[1], 'z')
                    + (a[1]/self.r) * self.differentiate(b[1], 'phi')
                    + a[2] * self.differentiate(b[1], 'r'))
                matderiv_r = ((-1*a[1]*b[1]/self.r)
                    + a[0] * self.differentiate(b[2], 'z')
                    + (a[1]/self.r) * self.differentiate(b[2], 'phi')
                    + a[2] * self.differentiate(b[2], 'r'))
                return np.array([matderiv_z, matderiv_phi, matderiv_r])
            if coordinates == 'Cartesian' or coordinates == 'cartesian':
                raise Exception("use spherical lol cartesian is borked")
        if b.shape == self.array_size:
            if coordinates == 'Spherical' or coordinates == 'spherical':
                matderiv = a[2] * (self.differentiate(b, 'phi') / (self.r * np.sin(self.theta))) + a[1] * (self.differentiate(b, 'theta') / self.r) + a[2] * self.differentiate(b, 'r')
                return matderiv
            if coordinates == 'Cylindrical' or coordinates == 'cylindrical':
                matderiv = a[2] * (self.differentiate(b, 'z')) + a[1] * (self.differentiate(b, 'phi') / self.r) + a[2] * self.differentiate(b, 'r')
                return matderiv
            if coordinates == 'Cartesian' or coordinates == 'cartesian':
                raise Exception("use spherical lol cartesian is borked")

    def radial_transport(self, q):
        """
        Computes radial transport of q calculated flux of material out of a cell from radial sides
        """
        self.get_primaries(get_vel_r=True)
        self.get_face_areas(get_r_face_areas=True)
        q_face = self.assign_to_face(q)
        vr_face = self.assign_to_face(self.vel_r)
        q_transport = np.zeros(self.array_size)
        q_transport = (q_face*vr_face*self.r_face_area)[:,:,:,1:] - (q_face*vr_face*self.r_face_area)[:,:,:,:-1]
        return q_transport

    def assign_to_face(self, q):
        """
        takes a cell center defined quantity and defines it on r faces by averaging between adjacent cells
        """
        self.native_grid()
        q_face = np.zeros(self.array_size + np.array([0,0,0,1]))
        q_face[:,:,:,1:-1] = (1/2)*(q[:,:,:,1:] + q[:,:,:,:-1])
        q_face[:,:,:,0] = q_face[:,:,:,1] - ((q_face[:,:,:,2]-q_face[:,:,:,1])/(self.r[:,:,:,2]-self.r[:,:,:,1])) * ((self.r[:,:,:,1]-self.r[:,:,:,0]))
        q_face[:,:,:,-1] = q_face[:,:,:,-2] + ((q_face[:,:,:,-2]-q_face[:,:,:,-3])/(self.r[:,:,:,-2]-self.r[:,:,:,-3])) * ((self.r[:,:,:,-1]-self.r[:,:,:,-2]))
        return q_face

    def alpha_visc(self, alpha):
        '''
        f = alpha (dP/dr phi_hat + 1/r dP/dphi r_hat)
        '''
        self.native_grid()
        self.get_primaries(get_press=True)

        scaling_factor = -1*(3/2)*alpha

        f_r = self.differentiate(self.press, 'phi') / self.r
        f_phi = self.differentiate(self.press, 'r') + (2*self.press / self.r)
        if self.gridtype == 'Spherical': #unconvinced this is true in spherical
            warnings.warn("check this is correct for spherical")
            f_unscaled = np.array([f_phi, np.zeros(self.array_size), f_r])
            return (scaling_factor* f_unscaled)
        if self.gridtype == 'Cylindrical':
            f_unscaled = np.array([np.zeros(self.array_size), f_phi, f_r])
            return (scaling_factor * f_unscaled)

    def alpha_torque(self, alpha):
        f = self.alpha_visc(alpha)
        torque_density = self.r*f[1] #z hat, there is only a z component to torque density
        return torque_density

    def get_boundary_flux(self, q, intermediates=False):
        """
        Takes the flux of a quantity out of the boundaries. Bassically q * vr * da integrated through theta and phi at the boundary.

        Parameters
        ----------
        q : list
            The quantity of which flux will be calculated
        intermediates : bool, default False
            Determines if intermediate values are returned. If false only the total boundary flux is returned if true, total, inner and outer boundary flux are returned in that order.
        """
        self.native_grid(get_r = True, get_phi = True, get_z = True)
        #flux at boundaries
        outer_flux = 0
        inner_flux = 0
        self.get_primaries(get_vel_r=True)
        '''
        inner_boundary_bool = (self.r_primitive[:].any() == self.possible_r[0])
        outer_boundary_bool = (self.r_primitive[:].any() == self.possible_r[-1])
        if self.gridtype == 'Spherical':
            in_boundary_thetaf = self.thetaf_primitive[inner_boundary_bool]
            in_boumdary_phif = self.phif_primitive[inner_boundary_bool]
            in_da = -(np.cos(in_boundary_thetaf[:,:+1]) - np.cos(in_boundary_theta[:,:])) * (in_boundary_phif[:,:+1] - in_boundary_phif[:,:]) * (self.possible_r[0] ** 2)
            inner_flux
        '''
        if self.gridtype == 'Spherical': 
            for n in range(self.NumMeshBlocks):
                for p in range(self.phi_len):
                    for t in range(self.theta_len):
                        for r in range(self.r_len):
                            if self.r_primitive[n,r] == self.possible_r[0]:
                                da = (np.cos(self.thetaf_primitive[n,t+1]) - np.cos(self.thetaf_primitive[n,t])) * (self.phif_primitive[n,p+1]-self.phif_primitive[n,p]) * (self.possible_r[0] ** 2)
                                inner_flux += -1 * (self.vel_r[n,p,t,r] * q[n,p,t,r] * da)
                            if self.r_primitive[n,r] == self.possible_r[-1]:
                                da = (np.cos(self.thetaf_primitive[n,t+1]) - np.cos(self.thetaf_primitive[n,t])) * (self.phif_primitive[n,p+1]-self.phif_primitive[n,p]) * (self.possible_r[-1] ** 2)
                                outer_flux += self.vel_r[n,p,t,r] * q[n,p,t,r] * da
        if self.gridtype == 'Cylindrical': 
            for n in range(self.NumMeshBlocks):
                for z in range(self.z_len):
                    for p in range(self.phi_len):
                        for r in range(self.r_len):
                            if self.r_primitive[n,r] == self.possible_r[0]:
                                da = (self.zf_primitive[n,z+1] - self.zf_primitive[n,z]) * (self.phif_primitive[n,p+1]-self.phif_primitive[n,p]) * (self.possible_r[0])
                                inner_flux += -1 * (self.vel_r[n,z,p,r] * q[n,z,p,r] * da)
                            if self.r_primitive[n,r] == self.possible_r[-1]:
                                da = (self.zf_primitive[n,z+1] - self.zf_primitive[n,z]) * (self.phif_primitive[n,p+1]-self.phif_primitive[n,p]) * (self.possible_r[-1])
                                outer_flux += self.vel_r[n,z,p,r] * q[n,z,p,r] * da
        total_flux = inner_flux + outer_flux
        if intermediates == True:
            return total_flux, inner_flux, outer_flux
        else:
            return total_flux

    def get_slice_flux(self, q, r_slicepoint):
        """
        Takes the flux of a quantity out of the boundaries. Bassically q * vr * da integrated through theta and phi at the boundary.

        Parameters
        ----------
        q : list
            The quantity of which flux will be calculated
        intermediates : bool, default False
            Determines if intermediate values are returned. If false only the total boundary flux is returned if true, total, inner and outer boundary flux are returned in that order.
        """
        self.native_grid(get_r = True, get_phi = True, get_z = True)
        #flux at boundaries
        flux = 0
        self.get_primaries(get_vel_r=True)
        '''
        inner_boundary_bool = (self.r_primitive[:].any() == self.possible_r[0])
        outer_boundary_bool = (self.r_primitive[:].any() == self.possible_r[-1])
        if self.gridtype == 'Spherical':
            in_boundary_thetaf = self.thetaf_primitive[inner_boundary_bool]
            in_boumdary_phif = self.phif_primitive[inner_boundary_bool]
            in_da = -(np.cos(in_boundary_thetaf[:,:+1]) - np.cos(in_boundary_theta[:,:])) * (in_boundary_phif[:,:+1] - in_boundary_phif[:,:]) * (self.possible_r[0] ** 2)
            inner_flux
        '''

        r_slicepoint = self.possible_r[np.argmin(abs(self.possible_r - r_slicepoint))]

        if self.gridtype == 'Spherical': 
            for n in range(self.NumMeshBlocks):
                for p in range(self.phi_len):
                    for t in range(self.theta_len):
                        for r in range(self.r_len):
                            if self.r_primitive[n,r] == r_slicepoint:
                                da = (np.cos(self.thetaf_primitive[n,t+1]) - np.cos(self.thetaf_primitive[n,t])) * (self.phif_primitive[n,p+1]-self.phif_primitive[n,p]) * (r_slicepoint ** 2)
                                flux += (self.vel_r[n,p,t,r] * q[n,p,t,r] * da)
        if self.gridtype == 'Cylindrical': 
            for n in range(self.NumMeshBlocks):
                for z in range(self.z_len):
                    for p in range(self.phi_len):
                        for r in range(self.r_len):
                            if self.r_primitive[n,r] == r_slicepoint:
                                da = (self.zf_primitive[n,z+1] - self.zf_primitive[n,z]) * (self.phif_primitive[n,p+1]-self.phif_primitive[n,p]) * (r_slicepoint)
                                flux += (self.vel_r[n,z,p,r] * q[n,z,p,r] * da)
        
        return flux


    def get_boundary_int(self, q, intermediates=False):
        """
        Takes the flux of a quantity at the boundaries. Bassically q * da integrated through theta and phi at the boundary.

        Parameters
        ----------
        q : list
            The quantity of which flux will be calculated
        intermediates : bool, default False
            Determines if intermediate values are returned. If false only the total boundary flux is returned if true, total, inner and outer boundary flux are returned in that order.
        """
        #flux at boundaries
        outer_flux = 0
        inner_flux = 0
        self.get_primaries(get_vel_r=True)
        '''
        inner_boundary_bool = (self.r_primitive[:].any() == self.possible_r[0])
        outer_boundary_bool = (self.r_primitive[:].any() == self.possible_r[-1])
        if self.gridtype == 'Spherical':
            in_boundary_thetaf = self.thetaf_primitive[inner_boundary_bool]
            in_boumdary_phif = self.phif_primitive[inner_boundary_bool]
            in_da = -(np.cos(in_boundary_thetaf[:,:+1]) - np.cos(in_boundary_theta[:,:])) * (in_boundary_phif[:,:+1] - in_boundary_phif[:,:]) * (self.possible_r[0] ** 2)
            inner_flux
        '''
        if self.gridtype == 'Spherical': 
            for n in range(self.NumMeshBlocks):
                for p in range(self.phi_len):
                    for t in range(self.theta_len):
                        for r in range(self.r_len):
                            if self.r_primitive[n,r] == self.possible_r[0]:
                                da = -(np.cos(self.thetaf_primitive[n,t+1]) - np.cos(self.thetaf_primitive[n,t])) * (self.phif_primitive[n,p+1]-self.phif_primitive[n,p]) * (self.possible_r[0] ** 2)
                                inner_flux += -1 * (q[n,p,t,r] * da)
                            if self.r_primitive[n,r] == self.possible_r[-1]:
                                da = -(np.cos(self.thetaf_primitive[n,t+1]) - np.cos(self.thetaf_primitive[n,t])) * (self.phif_primitive[n,p+1]-self.phif_primitive[n,p]) * (self.possible_r[-1] ** 2)
                                outer_flux += q[n,p,t,r] * da
        if self.gridtype == 'Cylindrical': 
            for n in range(self.NumMeshBlocks):
                for z in range(self.z_len):
                    for p in range(self.phi_len):
                        for r in range(self.r_len):
                            if self.r_primitive[n,r] == self.possible_r[0]:
                                da = -(self.zf_primitive[n,z+1] - self.zf_primitive[n,z]) * (self.phif_primitive[n,p+1]-self.phif_primitive[n,p]) * (self.possible_r[0])
                                inner_flux += -1 * (q[n,z,p,r] * da)
                            if self.r_primitive[n,r] == self.possible_r[-1]:
                                da = -(self.zf_primitive[n,z+1] - self.zf_primitive[n,z]) * (self.phif_primitive[n,p+1]-self.phif_primitive[n,p]) * (self.possible_r[-1])
                                outer_flux += q[n,z,p,r] * da
        total_flux = inner_flux + outer_flux
        if intermediates == True:
            return total_flux, inner_flux, outer_flux
        else:
            return total_flux

    def build_vectors(self):
        """
        Builds r, v, and r_hat vectors
        """
        if self.r_vec is None or self.v_vec is None or self.r_hat is None:
            if self.gridtype == 'Spherical':
                self.get_primaries(get_vel_r=True, get_vel_phi=True, get_vel_theta=True)
                self.native_grid(get_r=True)
                # rotating frame -> nonrotating frame
                vp = self.vel_phi + (self.r*np.sin(self.theta) * sim.orbital_Omega)
                # set up vectors
                self.r_vec = np.array([np.zeros(self.array_size), np.zeros(self.array_size), self.r])
                self.r_hat = np.array([np.zeros(self.array_size), np.zeros(self.array_size), np.full(self.array_size, 1)])
                self.v_vec = np.array([vp, self.vel_theta, self.vel_r])
                self.binary_omega = np.array([np.zeros(self.array_size), - np.full(self.array_size, sim.orbital_Omega) * np.sin(self.theta), np.full(self.array_size, sim.orbital_Omega) * np.cos(self.theta)])
            if self.gridtype == 'Cylindrical':
                self.get_primaries(get_vel_r=True, get_vel_phi=True, get_vel_z=True)
                self.native_grid(get_r=True, get_z=True)
                # rotating frame -> nonrotating frame
                vp = self.vel_phi + (self.r * sim.orbital_Omega)
                # set up vectors
                if self.NO_Z_GRAVITY:
                    self.r_vec = np.array([np.zeros(self.array_size), np.zeros(self.array_size), self.r])
                    self.v_vec = np.array([np.zeros(self.array_size), vp, self.vel_r])
                    self.r_hat = np.array([np.zeros(self.array_size), np.zeros(self.array_size), np.full(self.array_size, 1)])
                    self.binary_omega = np.array([np.full(self.array_size, sim.orbital_Omega), np.zeros(self.array_size), np.zeros(self.array_size)])
                else:
                    print("wtf2")
                    self.r_vec = np.array([self.z, np.zeros(self.array_size), self.r])
                    self.v_vec = np.array([self.vel_z, vp, self.vel_r])
                    self.r_hat = (1 / np.sqrt(2)) * np.array([np.full(self.array_size, 1), np.zeros(self.array_size), np.full(self.array_size, 1)])
                    self.binary_omega = np.array([np.full(self.array_size, sim.orbital_Omega), np.zeros(self.array_size), np.zeros(self.array_size)])

    def get_lrl(self, components=False):
        '''
        calculate LRL then returns its magnitude and direction as an angle from the x axis
        e = (v x r x v) / k - rhat where k = G M1.
        If components is True will recturn magnitude, angle from x axis, r component, phi component, x component, and y component in that order

        Parameters
        ----------
        components : bool
            Determings if componenents are also returned or if only magnitude and orientation are returned
        '''
    
        self._axes()
        lrl = np.zeros(self.vector_array_size)
        eccent = np.zeros(self.array_size)

        self.build_vectors()
        
        # Calculation
        for n in range(self.NumMeshBlocks):
            rxv = vec.ortho_cross(self.r_vec[:,n], self.v_vec[:,n])
            vxrxv = vec.ortho_cross(self.v_vec[:,n], rxv)
            lrl[:,n] = (vxrxv/sim.gm1) - self.r_hat[:,n]
            eccent[n] = vec.get_magnitude(lrl[:,n])

        cart_lrl = self.native_to_cart(lrl)

        if components==True:
            return lrl, cart_lrl
        
        if components == False:
            pos_x_lrl_bool = (cart_lrl[0] >= 0)
            pos_x_lrl_y = np.multiply(pos_x_lrl_bool, cart_lrl[1])
            neg_x_lrl_bool = np.logical_not(pos_x_lrl_bool)
            neg_x_lrl_y = np.multiply(neg_x_lrl_bool, cart_lrl[1])
            lrl_orient = (np.arctan(pos_x_lrl_y / cart_lrl[0]) % (2*np.pi)) + np.multiply(((np.pi + np.arctan(neg_x_lrl_y / cart_lrl[0])) % (2*np.pi)), neg_x_lrl_bool)
            return eccent, lrl_orient
        

    def get_C_vec(self, acceleration):
        """C = (a x (r x v) + v x (r x a)) / k where k = G M1.

        Returns:
            The time derivative of the rescaled LRL vector due to an
            acceleration in native coordinates
        """
        C = np.zeros(self.vector_array_size)

        self.build_vectors()
        a = acceleration
        
        # Calculation
        for n in range(self.NumMeshBlocks):
            rxv = vec.ortho_cross(self.r_vec[:,n], self.v_vec[:,n])
            axrxv = vec.ortho_cross(a[:,n], rxv)
            rxa = vec.ortho_cross(self.r_vec[:,n], a[:,n])
            vxrxa = vec.ortho_cross(self.v_vec[:,n], rxa)
            C[:,n] = np.array((axrxv + vxrxa) / sim.gm1)
        
        return C

    def native_to_cart(self, vector, flat=False):
        '''
        Converts vectors from native grid basis to cartesian basis.
        '''
        
        if self.gridtype == "Spherical":
            [p, t, r] = vector
            if flat == True:
                x = r*np.cos(self.phi) - p*np.sin(self.phi)
                y = r*np.sin(self.phi) + p*np.cos(self.phi)
                return np.array([x, y, 0])
            elif flat == False:
                x = r*np.sin(self.theta)*np.cos(self.phi) + t*np.cos(self.theta)*np.cos(self.phi) - p*np.sin(self.phi)
                y = r*np.sin(self.theta)*np.sin(self.phi) + t*np.cos(self.theta)*np.sin(self.phi) + p*np.cos(self.phi)
                z = r*np.cos(self.theta) - t*np.sin(self.theta)
                return np.array([x, y, z])
        if self.gridtype == "Cylindrical":
            self.native_grid(get_r=True, get_phi=True)
            [z, p, r] = vector
            x = r*np.cos(self.phi) - p*np.sin(self.phi)
            y = r*np.sin(self.phi) + p*np.cos(self.phi)
            return np.array([x, y, z])

#test = Athena_Analysis('/home/morgan/mnt/kitp/data/cvdisk/CVThick3/disk.out1.00640.athdf') #shapes (4720, 32, 4, 32)
#cv2 = Athena_Analysis('/home/morgan/mnt/kitp/data2/cvdisk/CVThin2/Data/disk.out1.00640.athdf') # shapes (4408, 96, 4, 48)
#test = Athena_Analysis('/home/morgan/mnt/kitp/data/cvdisk/superhump_3d_alpha03/disk.out1.02000.athdf')
#cyl1 = Athena_Analysis(file.data_loc+'Cyl_1/disk.out1.00002.athdf', grid_type="Cylindrical")