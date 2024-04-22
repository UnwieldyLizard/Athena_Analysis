from .roots.athena_analysis import *
from .roots.misc_func import *
import pickle

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

    def profile(self, fnum, plot=True, flatten=False, reynolds=False, just_torques=False, spread=True, vbound=1e2):
        self.torques = {}
        self.total_torques = {}
        if self.is_MHD:
            self.torque_list = ["tidal", "maxwell"]
        else:
            self.torque_list = ["tidal", "alpha"]
        if reynolds:
            self.torque_list.append("reynolds")

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
        reyn_rp = aa.rho*aa.vel_r*aa.vel_phi
        self.torques["reynold"] = -1 * aa.r * (aa.differentiate(reyn_rp, 'r')) - 2*(reyn_rp)
        if self.is_MHD:
            max_rp = -1*aa.B_r*aa.B_phi
            self.torques["maxwell"] = -1 * aa.r * (aa.differentiate(max_rp, 'r')) - 2*(max_rp)
        else:
            alpha_rp = (3/2)*self.alpha*aa.press
            self.torques["alpha"] = -1 * aa.r * (aa.differentiate(alpha_rp, 'r')) - 2*(alpha_rp)

        if flatten:
            self.sname = "flat"
            #Shell ave
            normalization_weight = aa.integrate(1, "shell")
            for key in self.torque_list:
                self.total_torques[key], self.torques[key], azimuthal_integral = aa.integrate(self.torques[key], "All", intermediates=True) #azimuthal will be ignored
                self.torques[key] = self.torques[key] / normalization_weight
            self.L_z = aa.integrate(self.L_z, "shell") / normalization_weight

            self.r_axis = aa.possible_r

            if plot:

                if just_torques:
                    if spread:
                        #plot
                        self.sname = "spread"
                        vert_num = 3
                        horz_num = 2
                        gs = gridspec.GridSpec(vert_num, horz_num)
                        fig = plt.figure(figsize=(horz_num*3, vert_num*3), dpi=300)
                        
                        ax_torques = fig.add_subplot(gs[0:2,:])
                        ax_image = fig.add_subplot(gs[2, 1])
                        ax_rho = fig.add_subplot(gs[2, 0])

                        for k, key in enumerate(self.torque_list): 
                            if k == 0: c = 0 #skipping the ugly orange at c=1 LOL
                            else: c = k + 1
                            ax_torques.plot(self.r_axis, self.torques[key], f"C{c}-", label=key)
                            #ax_torques.axhline(self.total_torques[key], c=f"C{c}", linestyle="dashed")
                        ax_torques.set_ylabel(r"$\tau$")
                        ax_torques.set_ylim([-1600,500])
                        ax_torques.set_xlabel("r")
                        ax_torques.set_title("Torques")
                        ax_torques.legend()

                        ax_rho.plot(self.r_axis, aa.integrate(aa.rho, "shell")/aa.integrate(1, "shell"))
                        ax_rho.set_xlabel("r")
                        ax_rho.set_ylabel(r"$\rho$")
                        ax_rho.set_ylim([-5,100])
                        ax_rho.set_title("Density")

                        rotation = -1*sim.orbital_Omega * aa.time
                        aa.midplane_colorplot(aa.rho, ax_image, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
                        ax_image.set_title("Density")

                        plt.subplots_adjust(top=(1-0.01*(20/vert_num)))
                        orbit = (self.time / sim.binary_period)
                        fig.suptitle(f"{self.dname} Orbit: {orbit:.2f}")
                        plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
                        plt.close()
                    else:
                        #plot
                        vert_num = 1
                        horz_num = 2
                        gs = gridspec.GridSpec(vert_num, horz_num)
                        fig = plt.figure(figsize=(horz_num*3, vert_num*3), dpi=300)
                        
                        ax_torques = fig.add_subplot(gs[0,1])
                        ax_rho = fig.add_subplot(gs[0,0])

                        for k, key in enumerate(self.torque_list): 
                            if k == 0: c = 0 #skipping the ugly orange at c=1 LOL
                            else: c = k + 1
                            ax_torques.plot(self.r_axis, self.torques[key], f"C{c}-", label=key)
                            #ax_torques.axhline(self.total_torques[key], c=f"C{c}", linestyle="dashed")
                        ax_torques.set_ylabel(r"$\tau$")
                        ax_torques.set_ylim([-1600,500])
                        ax_torques.set_xlabel("r")
                        ax_torques.set_title("Torques")
                        ax_torques.legend()

                        ax_rho.plot(self.r_axis, aa.integrate(aa.rho, "shell")/aa.integrate(1, "shell"))
                        ax_rho.set_xlabel("r")
                        ax_rho.set_ylabel(r"$\rho$")
                        ax_rho.set_ylim([-5,100])
                        ax_rho.set_title("Density")

                        plt.subplots_adjust(top=(1-0.01*(16/vert_num)))
                        orbit = (self.time / sim.binary_period)
                        fig.suptitle(f"{self.dname} Orbit: {orbit:.2f}")
                        plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
                        plt.close()
                else:
                    self.plot(fnum)

        else:
            self.sname = "full"
            self.total_torque = 0
            for key in self.torque_list:
                self.total_torque += self.torques[key]

            if plot:
                rotation = -1*sim.orbital_Omega * aa.time

                vert = 2
                horz = 3
                gs = gridspec.GridSpec(vert, horz)
                fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

                ax_rho = fig.add_subplot(gs[0, 0])
                ax_tid = fig.add_subplot(gs[1, 0])

                ax_alphamax = fig.add_subplot(gs[0, 1])
                ax_reyn = fig.add_subplot(gs[1, 1])

                ax_total = fig.add_subplot(gs[0, 2])

                aa.midplane_colorplot(aa.rho, ax_rho, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
                aa.midplane_colorplot(self.torques["tidal"], ax_tid, log=True, vbound=[-vbound,vbound], slicetype='z', rotation=rotation)
                ax_rho.set_title("Density")
                ax_tid.set_title("Tidal Torque")
                
                if self.is_MHD:
                    aa.midplane_colorplot(self.torques["maxwell"], ax_alphamax, log=True, vbound=[-vbound,vbound], slicetype='z', rotation=rotation)
                    ax_alphamax.set_title("Magnetic Torque")
                else:
                    aa.midplane_colorplot(self.torques["alpha"], ax_alphamax, log=True, vbound=[-vbound,vbound], slicetype='z', rotation=rotation)
                    ax_alphamax.set_title(r"$\alpha$ Torque")
                aa.midplane_colorplot(self.torques["reynold"], ax_reyn, log=True, vbound=[-vbound,vbound], slicetype='z', rotation=rotation)
                ax_reyn.set_title("Reynolds Torque")

                aa.midplane_colorplot(self.total_torque, ax_total, log=True, vbound=[-vbound,vbound], slicetype='z', rotation=rotation)
                ax_total.set_title("Total Torque")

                plt.tight_layout()
                orbit = (aa.time / sim.binary_period)
                plt.subplots_adjust(top=(1-0.01*(16/vert)))
                fig.suptitle(f"{self.dname} orbit: {orbit:.2f}")
                plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
                plt.close()

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
        self.sname = f"_averaged_{duration}"
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
        ax_torques.set_ylim([-100, 500])
        #ax_torques.set_ylim([-10000,10000])
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

    def alpha_analysis(self, start_fnum, duration):
        file_duration = duration * sim.filenums_per_orbit
        fnum_range = np.arange(start_fnum, start_fnum+file_duration, self.file_spacing)
        now = datetime.now()
        for i, fnum in enumerate(fnum_range):
            logging.info(f"fnum = {fnum}, {i}/{len(fnum_range)-1}")
            logging.info(datetime.now()-now)
            now = datetime.now()

            filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
            aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)
            aa.get_Bfields()
            Bpress = (aa.B_r*aa.B_r + aa.B_phi*aa.B_phi + aa.B_z*aa.B_z)/2
            aa.get_primaries(get_press=True)

            if i == 0:
                normalization_weight = aa.integrate(1, "shell")
                self.averaged_Maxwell = aa.integrate((-2/3)*(aa.B_r*aa.B_phi)/(aa.press+Bpress), "shell") / (normalization_weight*len(fnum_range))
                start_time = aa.time
                self.r_axis = aa.possible_r
            else:
                self.averaged_Maxwell += aa.integrate((-2/3)*(aa.B_r*aa.B_phi)/(aa.press+Bpress), "shell") / (normalization_weight*len(fnum_range))
                end_time = aa.time

        #plot
        times = [start_time, end_time]
        self.sname = f"_averaged_Balpha_{duration}"

        vert_num = 1
        horz_num = 2
        gs = gridspec.GridSpec(vert_num, horz_num)
        fig = plt.figure(figsize=(horz_num*3, vert_num*3), dpi=300)
        
        ax = fig.add_subplot(gs[0, 0])
        ax_zoom = fig.add_subplot(gs[0, 1])

        ax.plot(self.r_axis, self.averaged_Maxwell, "C2-", label=r"$\alpha$")
        ax.set_ylabel(r"$\alpha$")
        ax.set_title(r"$\alpha$ Due to Maxwell Stess")

        ax_zoom.plot(self.r_axis, self.averaged_Maxwell, "C2-", label=r"$\alpha$")
        ax.set_title(r"$\alpha$ Zoomed Way in")
        ax.set_ylim([0,1])

        plt.subplots_adjust(top=(1-0.01*(16/vert_num)))
        orbits = (times / sim.binary_period)
        fig.suptitle(f"Averaged from Orbit {orbits[0]:.2f} to {orbits[1]:.2f}")
        plt.legend()
        plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, start_fnum, self.sname))
        plt.close()

    def time_evolution(self, cutoff_radius, fnum_range, plot_every=10, step=None):
        self.sname = "_time_evolution"
        if step is None:
            step = self.file_spacing
        fnum_range = np.arange(fnum_range[0], fnum_range[-1]+1, step)
        now = datetime.now()
        orbits = np.zeros(len(fnum_range))
        cutoff_integrated_torque_max = np.zeros(len(fnum_range))
        cutoff_integrated_torque_tid = np.zeros(len(fnum_range))
        square_integrated_torque_max = np.zeros(len(fnum_range))
        square_integrated_torque_tid = np.zeros(len(fnum_range))
        alpha = np.zeros(len(fnum_range))
        for i, fnum in enumerate(fnum_range):
            logging.info(f"fnum = {fnum}")
            logging.info(datetime.now()-now)
            now = datetime.now()

            if i == 0 and fnum != 0:
                logging.info("setting up")
                found_load_point = False
                load_point = fnum - self.file_spacing
                while found_load_point == False:
                    if os.path.exists("%s/pickles/%s_pickle_%05d%s.dat" % (self.savedir, self.dname, load_point, self.sname)):
                        logging.info("Found data, loading from: %s" % load_point)
                        found_load_point = True
                    else:
                        load_point -= self.file_spacing
                    if load_point <= 0:
                        raise("No load point, you need to restart")
                with open("%s/pickles/%s_pickle_%05d%s.dat" % (self.savedir, self.dname, load_point, self.sname), "rb") as pickle_file:
                    data = pickle.load(pickle_file)
                i_0 = len(data["orbits"])
                orbits = np.zeros(len(np.arange(0, fnum_range[-1]+1, step)))
                cutoff_integrated_torque_max = np.zeros(len(np.arange(0, fnum_range[-1]+1, step)))
                cutoff_integrated_torque_tid = np.zeros(len(np.arange(0, fnum_range[-1]+1, step)))
                square_integrated_torque_max = np.zeros(len(np.arange(0, fnum_range[-1]+1, step)))
                square_integrated_torque_tid = np.zeros(len(np.arange(0, fnum_range[-1]+1, step)))
                alpha = np.zeros(len(np.arange(0, fnum_range[-1]+1, step)))
                orbits[:i_0] = data["orbits"]
                alpha[:i_0] = data["alpha"]
                cutoff_integrated_torque_max[:i_0] = data["cutoff_integrated_torque_max"]
                cutoff_integrated_torque_tid[:i_0] = data["cutoff_integrated_torque_tid"]
                square_integrated_torque_max[:i_0] = data["square_integrated_torque_max"]
                square_integrated_torque_tid[:i_0] = data["square_integrated_torque_tid"]
            elif i==0:
                i_0 = 0
            j = i + i_0
            
            filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
            aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)
            aa.get_Bfields()
            Bpress = (aa.B_r*aa.B_r + aa.B_phi*aa.B_phi + aa.B_z*aa.B_z)/2
            aa.get_primaries(get_press=True, get_rho=True)
            aa.get_potentials(get_companion_grav=True, get_accel=True)
            cutoff_idx = np.argmin(abs(aa.possible_r - cutoff_radius))

            tidal_torque = aa.r * (-1 * aa.rho * aa.gradient(aa.accel_pot + aa.companion_grav_pot, coordinates=self.grid_type))[1]
            maxwell_torque = aa.differentiate(aa.r*aa.r*aa.B_r*aa.B_phi, 'r') / aa.r
            tidal_torque = aa.integrate(tidal_torque, "shell")
            maxwell_torque = aa.integrate(maxwell_torque, "shell")
            orbits[j] = aa.time / sim.binary_period
            normalization_weight = aa.integrate(1, "All")
            alpha[j] = aa.integrate((-2/3)*(aa.B_r*aa.B_phi)/(aa.press+Bpress), "All") / (normalization_weight)
            cutoff_integrated_torque_max[j] = np.sum((maxwell_torque*aa.possible_dr_primitive)[cutoff_idx:])
            cutoff_integrated_torque_tid[j] = np.sum((tidal_torque*aa.possible_dr_primitive)[cutoff_idx:])
            square_integrated_torque_max[j] = np.sum(maxwell_torque*maxwell_torque*aa.possible_dr_primitive)
            square_integrated_torque_tid[j] = np.sum(tidal_torque*tidal_torque*aa.possible_dr_primitive)

            if (i % plot_every == 0):
                #pickle
                logging.info("pickling up to %s" % fnum)
                mkdir_if_not_exist("%s/pickles" % (self.savedir))
                with open("%s/pickles/%s_pickle_%05d%s.dat" % (self.savedir, self.dname, fnum, self.sname), "wb") as pickle_file:
                    pickle.dump({
                        "orbits": orbits[:j+1],
                        "alpha": alpha[:j+1],
                        "cutoff_integrated_torque_max": cutoff_integrated_torque_max[:j+1],
                        "cutoff_integrated_torque_tid": cutoff_integrated_torque_tid[:j+1],
                        "square_integrated_torque_max": square_integrated_torque_max[:j+1],
                        "square_integrated_torque_tid": square_integrated_torque_tid[:j+1],
                    }, pickle_file)
                #plot
                vert_num = 2
                horz_num = 2
                gs = gridspec.GridSpec(vert_num, horz_num)
                fig = plt.figure(figsize=(horz_num*3, vert_num*3), dpi=300)
                
                ax_alpha = fig.add_subplot(gs[0, 0])
                ax_cutoff = fig.add_subplot(gs[0, 1])
                ax_square = fig.add_subplot(gs[1, 1])

                ax_alpha.plot(orbits[:j+1], alpha[:j+1])
                ax_alpha.set_title(r"$\alpha$")
                ax_cutoff.plot(orbits[:j+1], cutoff_integrated_torque_tid[:j+1], label="Tid")
                ax_cutoff.plot(orbits[:j+1], cutoff_integrated_torque_max[:j+1], label="Max")
                ax_cutoff.set_title(f"total torque on disk outside r={cutoff_radius}")
                ax_cutoff.legend()
                ax_square.plot(orbits[:j+1], square_integrated_torque_tid[:j+1], label="Tid")
                ax_square.plot(orbits[:j+1], square_integrated_torque_max[:j+1], label="Max")
                ax_square.set_title("square integrated torques")
                ax_square.legend()

                plt.tight_layout()
                plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
                plt.close()


def compare_accretion(dnames, fnum_range, restart = False, plot_every=10, pickle_every=10):
    aname = "accretion"
    file_spacings = np.zeros(len(dnames), dtype=int)
    for d in range(len(dnames)):
        file_spacings[d] = find_file_spacing(dnames[d])
    savedir = file.savedir + "comparison" + "/" + aname
    pickldir = savedir + "/pickles"
    mkdir_if_not_exist(savedir)
    mkdir_if_not_exist(pickldir)
    file_spacing = np.gcd.reduce(file_spacings)

    if restart == False:
        logging.info("setting up")
        found_load_point = False
        load_point = fnum_range[0] - file_spacing
        while found_load_point == False:
            if os.path.exists("%s/pickles/pickle_%05d.dat" % (savedir, load_point)):
                logging.info("Found data, loading from: %s" % load_point)
                found_load_point = True
            else:
                load_point -= file_spacing
            if load_point <= 0:
                load_point = 0
                logging.info("No load point found, restarting")
                break
        if load_point != 0:
            with open("%s/pickles/pickle_%05d.dat" % (savedir, load_point), "rb") as pickle_file:
                data = pickle.load(pickle_file)
            offset = len(data["orbits"])
        else:
            offset = 0
    else:
        offset = 0

    orbits = np.zeros([(fnum_range[1] - fnum_range[0] + offset)])
    masses = np.zeros([len(dnames), (fnum_range[1] - fnum_range[0] + offset)])
    accretion_rates = np.zeros([len(dnames), 2 ,(fnum_range[1] - fnum_range[0] + offset)])
    #0 for inner boundary 1 for outer
    final_points = np.zeros(len(dnames))

    if (restart == False) and (dnames != data["dnames"]):
        raise("this functionality doesn't exist yet, please use consistent dnames")
    if (restart == False) and (load_point != 0):
        orbits[:offset] = data["orbits"]
        masses[:, :offset] = data["masses"][:,:offset]
        accretion_rates[:, : ,:offset] = data["accretion_rates"][:, : ,:offset]
        final_points = data["final_points"]

    for j, f in enumerate(np.arange(fnum_range[0], fnum_range[1], file_spacing)):
        i = j + offset
        vert = 1
        horz = 3
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
        
        ax_mass = fig.add_subplot(gs[0, 0])
        ax_mass.set_xlabel("orbits")
        #ax_mass.set_yscale("log")
        ax_mass.set_ylabel(r"M")
        ax_mass.set_title(r"Mass of Disk")

        ax_accr = fig.add_subplot(gs[0, 1])
        ax_accr.set_xlabel("orbits")
        #ax_accr.set_yscale("log")
        ax_accr.set_ylabel(r"$\dot{M}$")
        ax_accr.set_title(r"Accretion Rate")

        ax_flux = fig.add_subplot(gs[0, 2])
        ax_flux.set_xlabel("orbits")
        #ax_flux.set_yscale("log")
        ax_flux.set_ylabel(r"$\dot{M}$")
        ax_flux.set_title(r"Outer Boundary Flux")
        
        for d in range(len(dnames)):
            filename = "%s/disk.out1.%05d.athdf" % (file.data_loc + dnames[d], f)
            try:
                aa = Athena_Analysis(filename=filename, grid_type=file.grid_types[dnames[d]])
            except:
                if final_points[d] == 0:
                    logging.info(f"{dnames[d]} ended at f = {f}")
                    print(f"{dnames[d]} ended at f = {f}")
                    final_points[d] = i
                if np.all(final_points):
                    logging.info("All series terminated, pickling")
                    with open("%s/pickles/pickle_%05d.dat" % (savedir, f-1), "wb") as pickle_file:
                        pickle.dump({
                            "dnames": dnames,
                            "orbits": orbits[:i],
                            "masses": masses[:,:i],
                            "accretion_rates": accretion_rates[:,:,:i],
                            "final_points": final_points,
                        }, pickle_file)
                    return
            
            if (final_points[d] == 0):
                aa.get_primaries(get_rho=True)
                [total, inner, outer] = aa.get_boundary_flux(aa.rho, intermediates=True)

                masses[d, i] = aa.integrate(aa.rho, "All")
                accretion_rates[d, 0, i] = inner
                accretion_rates[d, 1, i] = outer

                orbits[i] = aa.time / sim.binary_period

                if i % plot_every == 0:
                    ax_mass.plot(orbits[:i+1], masses[d, :i+1], f"C{d}-", label=dnames[d])
                    ax_accr.plot(orbits[:i+1], accretion_rates[d, 0, :i+1], f"C{d}-", label=dnames[d])
                    ax_flux.plot(orbits[:i+1], accretion_rates[d, 1, :i+1], f"C{d}-", label=dnames[d])

            if final_points[d] != 0 and (i % plot_every == 0):
                ax_mass.plot(orbits[:int(final_points[d])], masses[d, :int(final_points[d])], f"C{d}-", label=dnames[d])
                ax_accr.plot(orbits[:int(final_points[d])], accretion_rates[d, 0, :int(final_points[d])], f"C{d}-", label=dnames[d])
                ax_flux.plot(orbits[:int(final_points[d])], accretion_rates[d, 1, :int(final_points[d])], f"C{d}-", label=dnames[d])

        if i % pickle_every == 0:
            with open("%s/pickles/pickle_%05d.dat" % (savedir, f), "wb") as pickle_file:
                pickle.dump({
                    "dnames": dnames,
                    "orbits": orbits[:i+1],
                    "masses": masses[:,:i+1],
                    "accretion_rates": accretion_rates[:,:,:i+1],
                    "final_points": final_points
                }, pickle_file)

        if i % plot_every == 0:
            plt.legend()
            plt.tight_layout()
            plt.savefig("%s/%s%05d.png" % (savedir, aname, f))
            ax_mass.set_yscale("log")
            ax_accr.set_yscale("log")
            ax_flux.set_yscale("log")
            plt.savefig("%s/%s_log%05d.png" % (savedir, aname, f))
            plt.close()

def append_dataset_accretion(new_dnames, fnum_range, plot_every=10, pickle_every = 100):
    aname = "accretion"
    file_spacings = np.zeros(len(new_dnames), dtype=int)
    for d in range(len(new_dnames)):
        file_spacings[d] = find_file_spacing(new_dnames[d])
    savedir = file.savedir + "comparison" + "/" + aname
    pickldir = savedir + "/pickles"
    mkdir_if_not_exist(savedir)
    mkdir_if_not_exist(pickldir)
    file_spacing = np.gcd.reduce(file_spacings)

    logging.info("setting up")
    found_load_point = False
    load_point = 10000 - file_spacing
    while found_load_point == False:
        if os.path.exists("%s/pickles/pickle_%05d.dat" % (savedir, load_point)):
            logging.info("Found data, loading from: %s" % load_point)
            found_load_point = True
        else:
            load_point -= file_spacing
        if load_point <= 0:
            load_point = 0
            logging.info("No load point found, restarting")
            break
    if load_point != 0:
        with open("%s/pickles/pickle_%05d.dat" % (savedir, load_point), "rb") as pickle_file:
            data = pickle.load(pickle_file)
        offset = len(data["orbits"])
    else:
        offset = 0

    old_dnames = data["dnames"]
    if offset > fnum_range[1]:
        max_orbit = offset
    else:
        max_orbit = fnum_range[1]
    dnames = np.concatenate([old_dnames, new_dnames])
    orbits = np.zeros([max_orbit])
    masses = np.zeros([len(dnames), (fnum_range[1] - fnum_range[0] + offset)])
    accretion_rates = np.zeros([len(dnames), 2 ,(fnum_range[1] - fnum_range[0] + offset)])
    final_points = np.zeros(len(dnames))
    
    orbits[:offset] = data["orbits"]
    masses[:len(old_dnames), :offset] = data["masses"][:, :offset]
    accretion_rates[:len(old_dnames), : ,:offset] = data["accretion_rates"][:, :, :offset]
    final_points[:len(old_dnames)] = data["final_points"]
    
    for j, f in enumerate(np.arange(fnum_range[0], fnum_range[1], file_spacing)):
        print(j, " file: ", f)
        i = j + offset
        vert = 1
        horz = 3
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
        
        ax_mass = fig.add_subplot(gs[0, 0])
        ax_mass.set_xlabel("orbits")
        #ax_mass.set_yscale("log")
        ax_mass.set_ylabel(r"M")
        ax_mass.set_title(r"Mass of Disk")

        ax_accr = fig.add_subplot(gs[0, 1])
        ax_accr.set_xlabel("orbits")
        #ax_accr.set_yscale("log")
        ax_accr.set_ylabel(r"$\dot{M}$")
        ax_accr.set_title(r"Accretion Rate")

        ax_flux = fig.add_subplot(gs[0, 2])
        ax_flux.set_xlabel("orbits")
        #ax_flux.set_yscale("log")
        ax_flux.set_ylabel(r"$\dot{M}$")
        ax_flux.set_title(r"Outer Boundary Flux")
        
        for d in range(len(dnames)):
            if d >= len(old_dnames):
                filename = "%s/disk.out1.%05d.athdf" % (file.data_loc + dnames[d], f)
                try:
                    aa = Athena_Analysis(filename=filename, grid_type=file.grid_types[dnames[d]])
                except:
                    if final_points[d] == 0:
                        logging.info(f"{dnames[d]} ended at f = {f}")
                        print(f"{dnames[d]} ended at f = {f}")
                        final_points[d] = i
                    if np.all(final_points[len(old_dnames):]):
                        logging.info("All new series terminated, pickling")
                        with open("%s/pickles/pickle_%05d.dat" % (savedir, max_orbit), "wb") as pickle_file:
                            pickle.dump({
                               "dnames": dnames,
                                "orbits": orbits[:i],
                                "masses": masses[:,:i],
                                "accretion_rates": accretion_rates[:,:,:i],
                                "final_points": final_points,
                            }, pickle_file)
                        return

                if (final_points[d] == 0):
                    aa.get_primaries(get_rho=True)
                    [total, inner, outer] = aa.get_boundary_flux(aa.rho, intermediates=True)

                    masses[d, i] = aa.integrate(aa.rho, "All")
                    accretion_rates[d, 0, i] = inner
                    accretion_rates[d, 1, i] = outer

                    orbits[i] = aa.time / sim.binary_period

                    if i % plot_every == 0:
                        ax_mass.plot(orbits[:i+1], masses[d, :i+1], f"C{d}-", label=dnames[d])
                        ax_accr.plot(orbits[:i+1], accretion_rates[d, 0, :i+1], f"C{d}-", label=dnames[d])
                        ax_flux.plot(orbits[:i+1], accretion_rates[d, 1, :i+1], f"C{d}-", label=dnames[d])

            if final_points[d] != 0 and (i % plot_every == 0):
                ax_mass.plot(orbits[:int(final_points[d])], masses[d, :int(final_points[d])], f"C{d}-", label=dnames[d])
                ax_accr.plot(orbits[:int(final_points[d])], accretion_rates[d, 0, :int(final_points[d])], f"C{d}-", label=dnames[d])
                ax_flux.plot(orbits[:int(final_points[d])], accretion_rates[d, 1, :int(final_points[d])], f"C{d}-", label=dnames[d])


            if (f % pickle_every == 0):
                with open("%s/pickles/pickle_%05d_new.dat" % (savedir, f), "wb") as pickle_file:
                    pickle.dump({
                        "dnames": dnames,
                        "orbits": orbits[:max_orbit],
                        "masses": masses[:,:max_orbit],
                        "accretion_rates": accretion_rates[:,:,:max_orbit],
                        "final_points": final_points
                    }, pickle_file)

        if i % plot_every == 0:
            plt.legend()
            plt.tight_layout()
            plt.savefig("%s/%s%05d.png" % (savedir, aname, f))
            ax_mass.set_yscale("log")
            ax_accr.set_yscale("log")
            ax_flux.set_yscale("log")
            plt.savefig("%s/%s_log%05d.png" % (savedir, aname, f))
            plt.close()