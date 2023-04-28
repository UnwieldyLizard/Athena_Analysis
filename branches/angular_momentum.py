from .roots.athena_analysis import *
from .roots.misc_func import *

logging.basicConfig(filename=file.logs_loc+"/angular_momentum.log", encoding='utf-8', level=logging.INFO)

class AngularMomentum():
    def __init__(self, dname):
        self.dname = dname
        self.aname = "_angular_momentum" #a for analysis
        self.sname = ""
        self.data_location = file.data_loc + dname
        self.grid_type = file.grid_types[dname]
        self.is_MHD = file.MHD[dname]
        self.alpha = file.alpha[dname]
        self.savedir = file.savedir + dname + "/" + dname + self.aname
        mkdir_if_not_exist(self.savedir)
        self.file_spacing = find_file_spacing(self.dname)
        
        self.time = None

    def profile(self, fnum, plot=True):
        self.torques = {}
        self.total_torques = {}
        if self.is_MHD:
            self.torque_list = ["tidal", "reynold", "maxwell"]
        else:
            self.torque_list = ["tidal", "reynold", "alpha"]
        
        filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
        aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)
        if self.time is not None:
            self.previous_time = self.time
        self.time = aa.time

        aa.get_primaries(get_rho=True, get_press=True, get_vel_r=True, get_vel_phi=True)
        aa.get_potentials(get_companion_grav=True, get_accel=True)
        if self.is_MHD:
            aa.get_Bfields()

        #angular momentum
        self.L_z = aa.r*aa.rho*aa.vel_phi

        #computing torque densities for z direction, ie for each force density (r x f)_z = r*f_phi. Recall vectors in [z, phi, r] order
        self.torques["tidal"] = aa.r * (-1 * aa.rho * aa.gradient(aa.accel_pot + aa.companion_grav_pot, coordinates=self.grid_type))[1]
        self.torques["reynold"] = -1 * aa.differentiate(aa.r*aa.r*aa.rho*aa.vel_phi*aa.vel_r, 'r') / aa.r
        if self.is_MHD:
            self.torques["maxwell"] = aa.differentiate(aa.r*aa.r*aa.B_r*aa.B_phi, 'r') / aa.r
        else:
            self.torques["alpha"] = -1 * aa.differentiate(aa.r*aa.r*(3/2)*self.alpha*aa.press, 'r') / aa.r

        #Shell ave
        normalization_weight = aa.integrate(1, "shell")
        for key in self.torque_list:
            self.total_torques[key], self.torques[key], azimuthal_integral = aa.integrate(self.torques[key], "All", intermediates=True) #azimuthal will be ignored
            self.torques[key] = self.torques[key] #/ normalization_weight
        self.L_z = aa.integrate(self.L_z, "shell") #/ normalization_weight

        self.r_axis = aa.possible_r

        if plot:
            #plot
            self.plot(fnum)

    def history(self, fnum_range, step=None):
        if step is None:
            step = self.file_spacing
        fnum_range = np.arange(fnum_range[0], fnum_range[-1]+1, step)
        now = datetime.now()
        for i, fnum in enumerate(fnum_range):
            logging.info(f"fnum = {fnum}")
            logging.info(datetime.now()-now)
            now = datetime.now()
            self.profile(fnum, plot=False)
            if i == 0 or i == 1:
                self.int_L_z = np.copy(self.L_z)
                dL_zdt = np.zeros(self.L_z.shape)
                for key in self.torque_list:
                    dL_zdt += self.torques[key]
            else:
                #self. quantities refer to previous timestep
                dL_zdt = np.zeros(self.L_z.shape)
                for key in self.torque_list:
                    dL_zdt += self.torques[key]
                self.int_L_z += (self.time - self.previous_time) * (dL_zdt + self.dL_zdt)/2

            #save values for use in next timestep
            self.dL_zdt = dL_zdt

            #plot
            self.plot(fnum, plot_int_L_z=True)

    def time_average_profile(self, start_fnum, duration):
        file_duration = duration * sim.filenums_per_orbit
        fnum_range = np.arange(start_fnum, start_fnum+file_duration, self.file_spacing)
        self.averaged_profile = {}
        self.averaged_totals = {}
        now = datetime.now()
        for i, fnum in enumerate(fnum_range):
            logging.info(f"fnum = {fnum}, {i}/{len(fnum_range)-1}")
            logging.info(datetime.now()-now)
            now = datetime.now()
            self.profile(fnum, plot=False)
            if i == 0:
                for key in self.torque_list:
                    self.averaged_profile[key] = self.torques[key] / len(fnum_range)
                    self.averaged_totals[key] = self.total_torques[key] / len(fnum_range)
                self.averaged_L_z = self.L_z / len(fnum_range)
                time = self.time
            else:
                for key in self.torque_list:
                    self.averaged_profile[key] += self.torques[key] / len(fnum_range)
                    self.averaged_totals[key] += self.total_torques[key] / len(fnum_range)
                self.averaged_L_z += self.L_z / len(fnum_range)

        #reassignment for plotting
        for key in self.torque_list:
            self.torques[key] = self.averaged_profile[key]
            self.total_torques[key] = self.averaged_totals[key]
        self.L_z = self.averaged_L_z

        #plot
        self.sname = "_averaged"
        self.plot(start_fnum, times=[time, self.time])
    
    def plot(self, fnum, times=None ,plot_int_L_z=False):
        #plot
        vert_num = 2
        horz_num = 1
        gs = gridspec.GridSpec(vert_num, horz_num)
        fig = plt.figure(figsize=(horz_num*3, vert_num*3), dpi=300)
        
        ax_L_z = fig.add_subplot(gs[0, 0])
        ax_torques = fig.add_subplot(gs[1, 0])

        ax_L_z.plot(self.r_axis, self.L_z, "C2-", label="measured")
        if plot_int_L_z:
            ax_L_z.plot(self.r_axis, self.int_L_z, "C9--", label="integrated")
        ax_L_z.set_ylabel(r"$L_z$")
        ax_L_z.set_ylim([1e2,1e6])
        ax_L_z.set_yscale("log")
        ax_L_z.set_title("Angular Momentum (log scale)")
        for k, key in enumerate(self.torque_list): 
            if k == 0: c = 0 #skipping the ugly orange at c=1 LOL
            else: c = k + 1
            ax_torques.plot(self.r_axis, self.torques[key], f"C{c}-", label=key+r" $\tau$="+f"{self.total_torques[key]:.2}")
            ax_torques.axhline(self.total_torques[key], c=f"C{c}", linestyle="dashed")
        ax_torques.set_ylabel(r"$\tau$")
        ax_torques.set_ylim([-10000,10000])
        ax_torques.set_xlabel("r")
        ax_torques.set_title("Torques")

        plt.subplots_adjust(top=(1-0.01*(16/vert_num)))
        if times is None:
            orbit = (self.time / sim.binary_period)
            fig.suptitle(f"Orbit: {orbit:.2f}")
        else:
            orbits = (times / sim.binary_period)
            fig.suptitle(f"Averaged from Orbit {orbits[0]:.2f} to {orbits[1]:.2f}")
        plt.legend()
        plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
        plt.close()
    