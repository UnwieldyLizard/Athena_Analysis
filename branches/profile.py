#Eventually needs to contain rewrite of mprofile code
from .roots.athena_analysis import *
from .roots.misc_func import *
import pickle


logging.basicConfig(filename=file.logs_loc+"/profile.log", encoding='utf-8', level=logging.INFO)

class Profile():
    def __init__(self, dname, start_point = 0):
        self.dname = dname
        self.aname = "_profile" #a for analysis
        self.sname = ""
        self.data_location = file.data_loc + dname
        self.grid_type = file.grid_types[dname]
        self.file_spacing = find_file_spacing(dname, start_point)
        self.is_MHD = file.MHD[dname]
        self.savedir = file.savedir + dname + "/" + dname + self.aname
        mkdir_if_not_exist(self.savedir)
        self.pickldir = self.savedir + "/pickles"
        mkdir_if_not_exist(self.pickldir)


    def temporal_profile(self, fnum_range, r_slicepoint = None, file_spacing=None, plot_every = 100):
        self.aname = "temporal"
        if file_spacing == None:
            file_spacing = self.file_spacing
        
        if r_slicepoint is not None:
            self.sname = "res"

        mass = np.zeros(fnum_range[1]-fnum_range[0])
        times = np.zeros(fnum_range[1]-fnum_range[0])

        for i, fnum in enumerate(range(fnum_range[0], fnum_range[1], file_spacing)):
            logging.info("fnum = %d" % fnum)

            filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)

            if not os.path.exists(filename):
                with open("%s/%s%s_pickle_%05d-%05d.dat" % (self.pickldir, self.dname, self.aname, fnum_range[0], fnum_range[1]), "wb") as pickle_file:
                    pickle.dump({"times": times, "mass": mass}, pickle_file)
                

            aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)

            aa.get_primaries(get_rho=True)
            aa.native_grid(get_r=True)

            if r_slicepoint is not None:
                r_idx = np.argmin(abs(aa.possible_r - r_slicepoint))

            times[i] = aa.time / sim.binary_period

            if r_slicepoint is None:
                mass[i] = aa.integrate(aa.rho, "all")
            else:
                mass[i] = aa.integrate(aa.rho, "phi", "Full", "z")[r_idx]

            if i % plot_every == 0:
                vert = 1
                horz = 1
                gs = gridspec.GridSpec(vert, horz)
                fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

                ax_rho = fig.add_subplot(gs[0, 0])
                ax_rho.plot(times[:i], mass[:i])
                
                if r_slicepoint is None:
                    ax_rho.set_title("Mass of Disk")
                    ax_rho.set_ylim([0, 130])
                else:
                    ax_rho.set_title(f"Mass of disk at r={aa.possible_r[r_idx]:.3}")
                    ax_rho.set_ylim([0, 5])
                ax_rho.set_xlabel("Binary Orbits")
                ax_rho.set_ylabel("M")

                plt.tight_layout()
                plt.subplots_adjust(top=(1-0.01*(16/vert)))
                fig.suptitle(self.dname)
                plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
                plt.close()


    def profile(self, fnum, r_slicepoint = 8):
        filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
        aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)

        aa.get_primaries(get_rho=True, get_vel_r=True, get_vel_phi=True, get_press=True)
        aa.get_face_areas(get_rcc_face_areas=True)
        aa.native_grid()
        if self.is_MHD:
            aa.get_Bfields()
        r_idx = np.argmin(abs(aa.possible_r - r_slicepoint))

        radial = {}
        vert = {}
        if self.is_MHD:
            keys = ["rho", "mdot", "rhovr", "press", "Bpress", "max_stress", "reyn_stress"]
        else:
            keys = ["rho", "mdot", "rhovr", "press", "reyn_stress"]

        if self.grid_type=="Cylindrical":
            vert_axis = aa.possible_z
        if self.grid_type=="Spherical":
            vert_axis = aa.possible_theta
        
        azavgvphi = aa.integrate(aa.vel_phi, "phi") / aa.integrate(1, "phi")
        turbulentvphi = np.zeros(aa.array_size) #defined as difference from mean phi velocity
        if self.grid_type == "Spherical":
            for n in range(aa.NumMeshBlocks):
                for k in range(aa.theta_len):
                    t = np.argwhere(aa.possible_theta == aa.theta_primitive[n, k])
                    for j in range(aa.r_len):
                        r = np.argwhere(aa.possible_r == aa.r_primitive[n, j])
                        turbulentvphi[n, :, k, j] = aa.vel_phi[n, :, k, j] - azavgvphi[t, r]
        if self.grid_type == "Cylindrical":
            for n in range(aa.NumMeshBlocks):
                for k in range(aa.z_len):
                    z = np.argwhere(aa.possible_z == aa.z_primitive[n, k])
                    for j in range(aa.r_len):
                        r = np.argwhere(aa.possible_r == aa.r_primitive[n, j])
                        turbulentvphi[n, k, :, j] = aa.vel_phi[n, k, :, j] - azavgvphi[z, r]

        radial["rho"], vert["rho"] = aa.integrate(aa.rho, "shell", intermediates=True)
        radial["mdot"], vert["mdot"] = aa.integrate(aa.rho * aa.vel_r * aa.rcc_face_area, "shell", intermediates=True)
        radial["rhovr"], vert["rhovr"] = aa.integrate(aa.rho * aa.vel_r, "shell", intermediates=True)
        radial["press"], vert["press"] = aa.integrate(aa.press, "shell", intermediates=True)
        radial["reyn_stress"], vert["reyn_stress"] = aa.integrate(aa.rho * aa.vel_r * turbulentvphi, variable="shell", intermediates=True)
        if self.is_MHD:
            radial["Bpress"], vert["Bpress"] = aa.integrate((aa.B_z**2+aa.B_phi**2+aa.B_z**2)/2, "shell", intermediates=True)
            radial["max_stress"], vert["max_stress"] = aa.integrate(-aa.B_r * aa.B_phi, "shell", intermediates=True)

        
        shell_normalization = aa.integrate(1, "shell")
        loop_normalization = (2*np.pi*aa.possible_r[r_idx])
        for key in keys:
            radial[key] = radial[key] / shell_normalization
            vert[key] = vert[key][:, r_idx] / loop_normalization

        radial["reyn_alpha"] = (2/3) * (radial["reyn_stress"] / radial["press"])
        vert["reyn_alpha"] = (2/3) * (vert["reyn_stress"] / vert["press"])
        if self.is_MHD:
            radial["max_alpha"] = (2/3) * (radial["max_stress"] / radial["press"])
            vert["max_alpha"] = (2/3) * (vert["max_stress"] / vert["press"])

        vert_num = 4
        horz_num = 2
        gs = gridspec.GridSpec(vert_num, horz_num)
        fig = plt.figure(figsize=(horz_num*3, vert_num*3), dpi=300)

        vertical_distribution = np.arange(0,40,1)
        
        ax_radrho = fig.add_subplot(gs[0, 0])
        ax_radpress = fig.add_subplot(gs[1, 0])
        ax_radalpha = fig.add_subplot(gs[2, 0])
        ax_radmdot = fig.add_subplot(gs[3, 0])

        ax_vertrho = fig.add_subplot(gs[0, 1])
        ax_vertpress = fig.add_subplot(gs[1, 1])
        ax_vertalpha = fig.add_subplot(gs[2, 1])
        ax_vertrhovr = fig.add_subplot(gs[3, 1])

        if self.grid_type=="Cylindrical":
            ax_radmdot.set_xlabel("r")
            ax_vertrhovr.set_xlabel("z, (r=%s)" % (r_slicepoint))
        if self.grid_type=="Spherical":
            ax_radmdot.set_xlabel("r")
            ax_vertrhovr.set_xlabel(r"$\theta$," + " (r=%s)" % (r_slicepoint))

        ax_radrho.plot(aa.possible_r, radial["rho"])
        ax_radmdot.plot(aa.possible_r, radial["mdot"])
        ax_radpress.plot(aa.possible_r, radial["press"], label="Gas")    
        ax_radalpha.plot(aa.possible_r, radial["reyn_alpha"], label="Reynolds")
        if self.is_MHD:
            ax_radpress.plot(aa.possible_r, radial["Bpress"], label="B")
            ax_radalpha.plot(aa.possible_r, radial["max_alpha"], label="Maxwell")

        ax_radrho.axvline(sim.three_one_res, c="red", linestyle="dashed", label="Resonant Radius")
        ax_radmdot.axvline(sim.three_one_res, c="red", linestyle="dashed", label="Resonant Radius")
        ax_radpress.axvline(sim.three_one_res, c="red", linestyle="dashed", label="Resonant Radius")
        ax_radalpha.axvline(sim.three_one_res, c="red", linestyle="dashed", label="Resonant Radius")

        ax_vertrho.plot(vert_axis, vert["rho"])
        ax_vertrhovr.plot(vert_axis, vert["rhovr"])
        ax_vertpress.plot(vert_axis, vert["press"], label="Gas")
        ax_vertalpha.plot(vert_axis, vert["reyn_alpha"], label="Reynolds")
        if self.is_MHD:
            ax_vertpress.plot(vert_axis, vert["Bpress"], label="B")
            ax_vertalpha.plot(vert_axis, vert["max_alpha"], label="Maxwell")

        ax_radrho.set_ylim([1e-3,1e3])
        ax_radrho.set_yscale("log")
        ax_radrho.legend()
        ax_radrho.set_title(r"Radial $\rho$ (log scale)")
        ax_vertrho.set_ylim([1e-3,1e3])
        ax_vertrho.set_title(r"Vertical $\rho$ (log scale)")
        ax_vertrho.set_yscale("log")
        ax_radmdot.set_title(r"Radial $\dot{M}$")
        ax_radmdot.set_ylim([-0.5,0.5])
        ax_vertrhovr.set_title(r"Vertical $\rho v_r$")
        ax_vertrhovr.set_ylim([-5,5])
        ax_radpress.set_title("Radial Pressures (log scale)")
        ax_radpress.set_yscale("log")
        ax_radpress.legend()
        ax_radpress.set_ylim([1e-2, 5e3])
        ax_vertpress.set_title(r"Vertical Pressures")
        ax_vertpress.legend()
        ax_vertpress.set_ylim([0,250])
        ax_radalpha.set_title(r"Radial $\alpha$ values (log scale)")
        ax_radalpha.set_yscale("log")
        ax_radalpha.legend()
        ax_radalpha.set_ylim([1e-6, 1e4])
        ax_vertalpha.set_title(r"Vertical $\alpha$ values (log scale)")
        ax_vertalpha.set_yscale("log")
        ax_vertalpha.legend()
        ax_vertalpha.set_ylim([1e-6, 1e4])

        plt.tight_layout()
        orbit = (aa.time / sim.binary_period)
        plt.subplots_adjust(top=(1-0.01*(16/vert_num)))
        fig.suptitle(f"orbit: {orbit:.2f}")
        plt.savefig("%s/%s%s%05d_r=%01d.png" % (self.savedir, self.dname, self.aname, fnum, r_slicepoint))
        plt.close()

    def plasma_beta(self, fnum):
        aname = "_plasma_beta" #a for analysis
        filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
        aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)

        aa.get_primaries(get_press=True)
        aa.get_Bfields()
        B_press = (aa.B_z**2 + aa.B_phi**2 + aa.B_r**2)/2

        plasma_beta = aa.press / B_press

        vert = 1
        horz = 2
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
        ax = fig.add_subplot(gs[0, 0])
        ax_vert = fig.add_subplot(gs[0, 1])

        aa.midplane_colorplot(plasma_beta, ax, vbound=[1e-4,1e4], log=True, cmap="seismic")
        aa.midplane_colorplot(plasma_beta, ax_vert, vbound=[1e-4,1e4], log=True, slicetype='y', cmap="seismic")

        orbit = (aa.time / sim.binary_period)
        plt.subplots_adjust(top=(1-0.01*(16/vert)))
        fig.suptitle(f"orbit: {orbit:.2f}")
        plt.savefig("%s/%s%s%05d.png" % (self.savedir, self.dname, aname, fnum))
        plt.close()

def compare_beta(dnames, fnum_range, inverse=False, restart = False, plot_every=10, pickle_every=10):
    if inverse:
        aname = "inverse_beta" #a for analysis
    else:
        aname = "beta"
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
    betas = np.zeros([len(dnames), (fnum_range[1] - fnum_range[0] + offset)])
    final_points = np.zeros(len(dnames))

    if (restart == False) and (dnames != data["dnames"]):
        raise("this functionality doesn't exist yet, please use consistent dnames")
    if (restart == False) and (load_point != 0):
        orbits[:offset] = data["orbits"]
        betas[:, :offset] = data["betas"][:,:offset]
        final_points = data["final_points"]

    for j, f in enumerate(np.arange(fnum_range[0], fnum_range[1], file_spacing)):
        i = j + offset
        vert = 1
        horz = 1
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
        ax = fig.add_subplot(gs[0, 0])
        ax.set_xlabel("orbits")
        ax.set_yscale("log")
        if inverse:
            ax.set_ylabel(r"$\frac{1}{\beta}$")
            ax.set_title(r"Mass Weighted Average $\frac{1}{\beta}$")
        else:
            ax.set_ylabel(r"$\beta$")
            ax.set_title(r"Mass Weighted Average $\beta$")
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
                            "betas": betas[:,:i],
                            "final_points": final_points,
                        }, pickle_file)
                    return
            
            if (final_points[d] == 0):
                aa.get_primaries(get_press=True, get_rho=True)
                aa.get_Bfields()
                B_press = (aa.B_z**2 + aa.B_phi**2 + aa.B_r**2)/2

                if inverse:
                    plasma_beta = B_press / aa.press
                else:
                    plasma_beta = aa.press / B_press
                mass = aa.integrate(aa.rho, "All")

                betas[d, i] = aa.integrate((aa.rho/mass) * plasma_beta, "All")
                orbits[i] = aa.time / sim.binary_period

                if i % plot_every == 0:
                    ax.plot(orbits[:i+1], betas[d, :i+1], f"C{d}-", label=dnames[d])

            if final_points[d] != 0 and (i % plot_every == 0):
                ax.plot(orbits[:int(final_points[d])], betas[d,:int(final_points[d])], f"C{d}-", label=dnames[d])

            if i % pickle_every == 0:
                with open("%s/pickles/pickle_%05d.dat" % (savedir, f), "wb") as pickle_file:
                    pickle.dump({
                        "inverse": inverse,
                        "dnames": dnames,
                        "orbits": orbits[:i+1],
                        "betas": betas[:,:i+1],
                        "final_points": final_points
                    }, pickle_file)

        if i % plot_every == 0:
            plt.legend()
            plt.savefig("%s/%s%05d.png" % (savedir, aname, f))
            plt.close()

def get_late_betas(dnames, cutoff_orbit, load_point, inverse = False):
    if inverse:
        aname = "inverse_beta" #a for analysis
    else:
        aname = "beta"
    file_spacings = np.zeros(len(dnames), dtype=int)
    for d in range(len(dnames)):
        file_spacings[d] = find_file_spacing(dnames[d])
    savedir = file.savedir + "comparison" + "/" + aname
    pickldir = savedir + "/pickles"
    mkdir_if_not_exist(savedir)
    mkdir_if_not_exist(pickldir)
    file_spacing = np.gcd.reduce(file_spacings)

    logging.info("setting up")
    found_load_point = False
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

    orbits = data["orbits"]
    betas = data["betas"]
    final_points = data["final_points"]

    vert = 1
    horz = 1
    gs = gridspec.GridSpec(vert, horz)
    fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlabel("orbits")
    ax.set_yscale("log")
    ax.set_ylim([10,1e4])
    if inverse:
        ax.set_ylabel(r"$\frac{1}{\beta}$")
        ax.set_title(r"Mass Weighted Average $\frac{1}{\beta}$")
    else:
        ax.set_ylabel(r"$\beta$")
        ax.set_title(r"Mass Weighted Average $\beta$")

    cutoff_idx = np.argmin(abs(orbits - cutoff_orbit))
    late_ave_betas = np.zeros(len(dnames))
    for d in range(len(dnames)):
        if final_points[d] != 0:
            late_time_nomalization = len(betas[d, cutoff_idx:int(final_points[d])])
            late_ave_betas[d] = sum(betas[d, cutoff_idx:int(final_points[d])]) / late_time_nomalization
            print(dnames[d] + ": " + str(late_ave_betas[d]))
            ax.plot(orbits[:int(final_points[d])], betas[d,:int(final_points[d])], f"C{d}-", label=dnames[d])
        else:
            late_time_nomalization = len(betas[d, cutoff_idx:])
            late_ave_betas[d] = sum(betas[d, cutoff_idx:]) / late_time_nomalization
            print(dnames[d] + ": " + str(late_ave_betas[d]))
            ax.plot(orbits[:], betas[d,:], f"C{d}-", label=dnames[d])

    plt.legend()
    plt.savefig("%s/%s%05d.png" % (savedir, aname, load_point))
    plt.close()



def compare_alpha(dnames, fnum_range, plot_every=10, pickle_every = 100, res = False, B_field_focus = False):
    if res:
        aname = "alpha_res"
    else:
        aname = "alpha"
    if B_field_focus:
        aname += "_B"
    file_spacings = np.zeros(len(dnames), dtype=int)
    for d in range(len(dnames)):
        file_spacings[d] = find_file_spacing(dnames[d])
    savedir = file.savedir + "comparison" + "/" + aname
    pickldir = savedir + "/pickles"
    mkdir_if_not_exist(savedir)
    mkdir_if_not_exist(pickldir)
    file_spacing = np.gcd.reduce(file_spacings)

    logging.info("setting up")
    found_load_point = False
    load_point = fnum_range[0]
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

    orbits = np.zeros([(fnum_range[1] - fnum_range[0] + offset)])
    alphas = np.zeros([len(dnames), (fnum_range[1] - fnum_range[0] + offset), 4])
    final_points = np.zeros(len(dnames))

    if load_point != 0:
        if dnames != data["dnames"]:
            raise("this functionality doesn't exist yet, please use consistent dnames")
        orbits[:offset] = data["orbits"]
        alphas[:, :offset] = data["alphas"][:,:offset]
        final_points = data["final_points"]

    def mass_weighted_ave(aa, value):
        aa.get_primaries(get_rho=True)
        if res:
            flat = aa.integrate(aa.rho * value, "z")
            flat_rho = aa.integrate(aa.rho, "z")
            res_r_idx = np.argmin(abs(aa.possible_r - sim.three_one_res))
            res_loop = flat[:, res_r_idx]
            res_loop_rho = flat_rho[:, res_r_idx]
            ave = np.sum(res_loop) / np.sum(res_loop_rho)
        else:
            ave = aa.integrate(aa.rho * value, "All") / aa.integrate(aa.rho, "All")
        return ave

    for j, f in enumerate(np.arange(fnum_range[0], fnum_range[1], file_spacing)):
        print(j, " file: ", f)
        i = j + offset
        if f % plot_every == 0:
            vert = 2
            horz = 2
            gs = gridspec.GridSpec(vert, horz)
            fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
            
            ax_max = fig.add_subplot(gs[0, 0])
            ax_max.set_xlabel("orbits")
            ax_max.set_ylabel(r"$\alpha$")
            ax_max.set_title(r"$\alpha_M$")
            ax_max.set_ylim([-1,1])

            ax_freyn = fig.add_subplot(gs[1, 0])
            ax_freyn.set_xlabel("orbits")
            if B_field_focus:
                ax_freyn.set_ylabel(r"P_B")
                ax_freyn.set_title("Magnetic Pressure")
            else:
                ax_freyn.set_ylabel(r"$\alpha$")
                ax_freyn.set_title(r"$\rho v_r v_{\phi}$")

            ax_treyn = fig.add_subplot(gs[0, 1])
            ax_treyn.set_xlabel("orbits")
            if B_field_focus:
                ax_treyn.set_ylabel(r"$B_z$")
                ax_treyn.set_title(r"$B_z$")
            else:
                ax_treyn.set_ylabel(r"$\alpha$")
                ax_treyn.set_title(r"$\rho \left(v_r - \left<v_r\right>\right) \left(v_{\phi} - \left<v_{\phi}\right>\right)$")

            if not B_field_focus:
                ax_creyn = fig.add_subplot(gs[1, 1])
                ax_creyn.set_xlabel("orbits")
                ax_creyn.set_ylabel(r"$\alpha$")
                ax_creyn.set_title(r"$\rho v_r \left(v_{\phi} - \left<v_{\phi}\right>\right)$")
        
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
                            "alphas": alphas[:,:i],
                            "final_points": final_points,
                        }, pickle_file)
                    return

            if final_points[d] == 0:
                aa.get_primaries(get_press=True, get_rho=True, get_vel_phi=True, get_vel_r=True)
                aa.get_Bfields()
                
                if B_field_focus:
                        maxwell_stress = (aa.B_r*aa.B_phi)
                        B_press = (aa.B_z**2 + aa.B_phi**2 + aa.B_r**2)/2
                        average_press = mass_weighted_ave(aa, aa.press)
                        alphas[d, i, 0] = (-2/3)*mass_weighted_ave(aa, maxwell_stress) / average_press
                        alphas[d, i, 1] = mass_weighted_ave(aa, aa.B_z)
                        alphas[d, i, 2] = mass_weighted_ave(aa, B_press)
                else:
                    azavgvphi = aa.integrate(aa.rho * aa.vel_phi, "phi") / aa.integrate(aa.rho, "phi")
                    turbulentvphi = np.zeros(aa.array_size)
                    azavgvr = aa.integrate(aa.rho * aa.vel_r, "phi") / aa.integrate(aa.rho, "phi")
                    turbulentvr = np.zeros(aa.array_size)
                    if file.grid_types[dnames[d]] == "Cylindrical":
                        for n in range(aa.NumMeshBlocks):
                            for k in range(aa.z_len):
                                z = np.argwhere(aa.possible_z == aa.z_primitive[n, k])
                                for j in range(aa.r_len):
                                    r = np.argwhere(aa.possible_r == aa.r_primitive[n, j])
                                    turbulentvphi[n, k, :, j] = aa.vel_phi[n, k, :, j] - azavgvphi[z, r]
                                    turbulentvr[n, k, :, j] = aa.vel_r[n, k, :, j] - azavgvr[z, r]
                    
                    maxwell_stress = (aa.B_r*aa.B_phi)
                    turbulent_stress = (aa.rho * turbulentvr * turbulentvphi)
                    full_reyn = (aa.rho*aa.vel_r*aa.vel_phi)
                    circl_reyn = (aa.rho*aa.vel_r*turbulentvphi)
                    average_press = mass_weighted_ave(aa, aa.press)

                    alphas[d, i, 0] = (-2/3)*mass_weighted_ave(aa, maxwell_stress) / average_press
                    alphas[d, i, 1] = (-2/3)*mass_weighted_ave(aa, turbulent_stress) / average_press
                    alphas[d, i, 2] = (-2/3)*mass_weighted_ave(aa, full_reyn) / average_press
                    alphas[d, i, 3] = (-2/3)*mass_weighted_ave(aa, circl_reyn) / average_press
                orbits[i] = aa.time / sim.binary_period

            if f % plot_every == 0:
                if final_points[d] == 0:
                    ax_max.plot(orbits[:i+1], alphas[d, :i+1, 0], f"C{d}-", label=dnames[d])
                    ax_treyn.plot(orbits[:i+1], alphas[d, :i+1, 1], f"C{d}-", label=dnames[d])
                    ax_freyn.plot(orbits[:i+1], alphas[d, :i+1, 2], f"C{d}-", label=dnames[d])
                    if not B_field_focus:
                        ax_creyn.plot(orbits[:i+1], alphas[d, :i+1, 3], f"C{d}-", label=dnames[d])
                else:
                    ax_max.plot(orbits[:int(final_points[d])], alphas[d, :int(final_points[d]), 0], f"C{d}-", label=dnames[d])
                    ax_treyn.plot(orbits[:int(final_points[d])], alphas[d, :int(final_points[d]), 1], f"C{d}-", label=dnames[d])
                    ax_freyn.plot(orbits[:int(final_points[d])], alphas[d, :int(final_points[d]), 2], f"C{d}-", label=dnames[d])
                    if not B_field_focus:
                        ax_creyn.plot(orbits[:int(final_points[d])], alphas[d, :int(final_points[d]), 3], f"C{d}-", label=dnames[d])

            if f % pickle_every == 0:
                with open("%s/pickles/pickle_%05d.dat" % (savedir, f), "wb") as pickle_file:
                    pickle.dump({
                        "dnames": dnames,
                        "orbits": orbits[:i+1],
                        "alphas": alphas[:,:i+1],
                        "final_points" : final_points,
                    }, pickle_file)

        if f % plot_every == 0:
            plt.legend()
            plt.subplots_adjust(top=(1-0.01*(16/vert)))
            if B_field_focus:
                title = "Mass Weighted Average quantities"
            else:
                title = r"Mass Weighted Average $\alpha$s"
            if res:
                title += " at Resonant Radius"
            plt.suptitle(title)
            plt.tight_layout()
            plt.savefig("%s/%s%05d.png" % (savedir, aname, f))
            plt.close()

def append_dataset_alpha(new_dnames, fnum_range, plot_every=10, pickle_every = 100, res = False, B_field_focus = False):
    if res:
        aname = "alpha_res"
    else:
        aname = "alpha"
    if B_field_focus:
        aname += "_B"
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
    load_point = 2500 - file_spacing
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
    alphas = np.zeros([len(dnames), max_orbit, 4])
    final_points = np.zeros(len(dnames))
    
    orbits[:offset] = data["orbits"]
    alphas[:len(old_dnames), :offset, :] = data["alphas"][:,:offset, :]
    final_points[:len(old_dnames)] = data["final_points"]

    def mass_weighted_ave(aa, value):
        aa.get_primaries(get_rho=True)
        if res:
            flat = aa.integrate(aa.rho * value, "z")
            flat_rho = aa.integrate(aa.rho, "z")
            res_r_idx = np.argmin(abs(aa.possible_r - sim.three_one_res))
            res_loop = flat[:, res_r_idx]
            res_loop_rho = flat_rho[:, res_r_idx]
            ave = np.sum(res_loop) / np.sum(res_loop_rho)
        else:
            ave = aa.integrate(aa.rho * value, "All") / aa.integrate(aa.rho, "All")
        return ave
    
    for j, f in enumerate(np.arange(fnum_range[0], fnum_range[1], file_spacing)):
        print(j, " file: ", f)
        i = j #+ offset
        if f % plot_every == 0:
            vert = 2
            horz = 2
            gs = gridspec.GridSpec(vert, horz)
            fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
            
            ax_max = fig.add_subplot(gs[0, 0])
            ax_max.set_xlabel("orbits")
            ax_max.set_ylabel(r"$\alpha$")
            ax_max.set_title(r"$\alpha_M$")

            ax_freyn = fig.add_subplot(gs[1, 0])
            ax_freyn.set_xlabel("orbits")
            if B_field_focus:
                ax_freyn.set_ylabel(r"P_B")
                ax_freyn.set_title("Magnetic Pressure")
            else:
                ax_freyn.set_ylabel(r"$\alpha$")
                ax_freyn.set_title(r"$\rho v_r v_{\phi}$")

            ax_treyn = fig.add_subplot(gs[0, 1])
            ax_treyn.set_xlabel("orbits")
            if B_field_focus:
                ax_treyn.set_ylabel(r"$B_z$")
                ax_treyn.set_title(r"$B_z$")
            else:
                ax_treyn.set_ylabel(r"$\alpha$")
                ax_treyn.set_title(r"$\rho \left(v_r - \left<v_r\right>\right) \left(v_{\phi} - \left<v_{\phi}\right>\right)$")

            if not B_field_focus:
                ax_creyn = fig.add_subplot(gs[1, 1])
                ax_creyn.set_xlabel("orbits")
                ax_creyn.set_ylabel(r"$\alpha$")
                ax_creyn.set_title(r"$\rho v_r \left(v_{\phi} - \left<v_{\phi}\right>\right)$")
        
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
                                "orbits": orbits[:max_orbit],
                                "alphas": alphas[:,:max_orbit],
                                "final_points": final_points,
                            }, pickle_file)
                        vert = 2
                        horz = 2
                        gs = gridspec.GridSpec(vert, horz)
                        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
                        
                        ax_max = fig.add_subplot(gs[0, 0])
                        ax_max.set_xlabel("orbits")
                        ax_max.set_ylabel(r"$\alpha$")
                        ax_max.set_title(r"$\alpha_M$")

                        ax_freyn = fig.add_subplot(gs[1, 0])
                        ax_freyn.set_xlabel("orbits")
                        if B_field_focus:
                            ax_freyn.set_ylabel(r"P_B")
                            ax_freyn.set_title("Magnetic Pressure")
                        else:
                            ax_freyn.set_ylabel(r"$\alpha$")
                            ax_freyn.set_title(r"$\rho v_r v_{\phi}$")

                        ax_treyn = fig.add_subplot(gs[0, 1])
                        ax_treyn.set_xlabel("orbits")
                        if B_field_focus:
                            ax_treyn.set_ylabel(r"$B_z$")
                            ax_treyn.set_title(r"$B_z$")
                        else:
                            ax_treyn.set_ylabel(r"$\alpha$")
                            ax_treyn.set_title(r"$\rho \left(v_r - \left<v_r\right>\right) \left(v_{\phi} - \left<v_{\phi}\right>\right)$")

                        ax_creyn = fig.add_subplot(gs[1, 1])
                        ax_creyn.set_xlabel("orbits")
                        ax_creyn.set_ylabel(r"$\alpha$")
                        ax_creyn.set_title(r"$\rho v_r \left(v_{\phi} - \left<v_{\phi}\right>\right)$")
                        for dname in range(len(dnames)):
                            ax_max.plot(int(final_points[dname]), alphas[dname, :int(final_points[dname]), 0], f"C{dname}-", label=dnames[dname])
                            ax_treyn.plot(int(final_points[dname]), alphas[dname, :int(final_points[dname]), 1], f"C{dname}-", label=dnames[dname])
                            ax_freyn.plot(int(final_points[dname]), alphas[dname, :int(final_points[dname]), 2], f"C{dname}-", label=dnames[dname])
                            ax_creyn.plot(int(final_points[dname]), alphas[dname, :int(final_points[dname]), 3], f"C{dname}-", label=dnames[dname])
                        plt.legend()
                        plt.subplots_adjust(top=(1-0.01*(16/vert)))
                        if res:
                            fig.suptitle(r"Mass Weighted Average $\alpha$s at Resonant Radius")
                        else:
                            fig.suptitle(r"Mass Weighted Average $\alpha$s")
                        plt.tight_layout()
                        plt.savefig("%s/%s%05d.png" % (savedir, aname, f))
                        plt.close()
                        return

                aa.get_primaries(get_press=True, get_rho=True, get_vel_phi=True, get_vel_r=True)
                aa.get_Bfields()

                if final_points[d] == 0:
                    if B_field_focus:
                        maxwell_stress = (aa.B_r*aa.B_phi)
                        B_press = (aa.B_z**2 + aa.B_phi**2 + aa.B_r**2)/2
                        average_press = mass_weighted_ave(aa, aa.press)
                        alphas[d, i, 0] = (-2/3)*mass_weighted_ave(aa, maxwell_stress) / average_press
                        alphas[d, i, 1] = mass_weighted_ave(aa, aa.B_z)
                        alphas[d, i, 2] = mass_weighted_ave(aa, B_press)
                    else:
                        azavgvphi = aa.integrate(aa.rho * aa.vel_phi, "phi") / aa.integrate(aa.rho, "phi")
                        turbulentvphi = np.zeros(aa.array_size)
                        azavgvr = aa.integrate(aa.rho * aa.vel_r, "phi") / aa.integrate(aa.rho, "phi")
                        turbulentvr = np.zeros(aa.array_size)
                        if file.grid_types[dnames[d]] == "Cylindrical":
                            for n in range(aa.NumMeshBlocks):
                                for k in range(aa.z_len):
                                    z = np.argwhere(aa.possible_z == aa.z_primitive[n, k])
                                    for j in range(aa.r_len):
                                        r = np.argwhere(aa.possible_r == aa.r_primitive[n, j])
                                        turbulentvphi[n, k, :, j] = aa.vel_phi[n, k, :, j] - azavgvphi[z, r]
                                        turbulentvr[n, k, :, j] = aa.vel_r[n, k, :, j] - azavgvr[z, r]
                        
                        maxwell_stress = (aa.B_r*aa.B_phi)
                        turbulent_stress = (aa.rho * turbulentvr * turbulentvphi)
                        full_reyn = (aa.rho*aa.vel_r*aa.vel_phi)
                        circl_reyn = (aa.rho*aa.vel_r*turbulentvphi)
                        average_press = mass_weighted_ave(aa, aa.press)

                        alphas[d, i, 0] = (-2/3)*mass_weighted_ave(aa, maxwell_stress) / average_press
                        alphas[d, i, 1] = (-2/3)*mass_weighted_ave(aa, turbulent_stress) / average_press
                        alphas[d, i, 2] = (-2/3)*mass_weighted_ave(aa, full_reyn) / average_press
                        alphas[d, i, 3] = (-2/3)*mass_weighted_ave(aa, circl_reyn) / average_press
                    orbits[i] = aa.time / sim.binary_period

            if f % plot_every == 0:
                if final_points[d] == 0:
                    ax_max.plot(orbits[:i+1], alphas[d, :i+1, 0], f"C{d}-", label=dnames[d])
                    ax_treyn.plot(orbits[:i+1], alphas[d, :i+1, 1], f"C{d}-", label=dnames[d])
                    ax_freyn.plot(orbits[:i+1], alphas[d, :i+1, 2], f"C{d}-", label=dnames[d])
                    if not B_field_focus:
                        ax_creyn.plot(orbits[:i+1], alphas[d, :i+1, 3], f"C{d}-", label=dnames[d])
                else:
                    ax_max.plot(orbits[:int(final_points[d])], alphas[d, :int(final_points[d]), 0], f"C{d}-", label=dnames[d])
                    ax_treyn.plot(orbits[:int(final_points[d])], alphas[d, :int(final_points[d]), 1], f"C{d}-", label=dnames[d])
                    ax_freyn.plot(orbits[:int(final_points[d])], alphas[d, :int(final_points[d]), 2], f"C{d}-", label=dnames[d])
                    if not B_field_focus:
                        ax_creyn.plot(orbits[:int(final_points[d])], alphas[d, :int(final_points[d]), 3], f"C{d}-", label=dnames[d])


            if (f % pickle_every == 0):
                with open("%s/pickles/pickle_%05d_new.dat" % (savedir, f), "wb") as pickle_file:
                    pickle.dump({
                        "dnames": dnames,
                        "orbits": orbits[:max_orbit],
                        "alphas": alphas[:,:max_orbit],
                        "final_points" : final_points,
                    }, pickle_file)

        if f % plot_every == 0:
            plt.legend()
            plt.subplots_adjust(top=(1-0.01*(16/vert)))
            if B_field_focus:
                title = "Mass Weighted Average quantities"
            else:
                title = r"Mass Weighted Average $\alpha$s"
            if res:
                title += " at Resonant Radius"
            fig.suptitle(title)
            plt.tight_layout()
            plt.savefig("%s/%s%05d.png" % (savedir, aname, f))
            plt.close()

def replot_compare_alpha(fnum, res = False, B_field_focus = False):
    if res:
        aname = "alpha_res"
    else:
        aname = "alpha"
    if B_field_focus:
        aname += "_B"
    savedir = file.savedir + "comparison" + "/" + aname
    with open("%s/pickles/pickle_%05d.dat" % (savedir, fnum), "rb") as pickle_file:
        data = pickle.load(pickle_file)

    dnames = data["dnames"]
    alphas = data["alphas"]
    orbits = data["orbits"]
    final_points = data["final_points"]

    vert = 2
    horz = 2
    gs = gridspec.GridSpec(vert, horz)
    fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
    
    ax_max = fig.add_subplot(gs[0, 0])
    ax_max.set_xlabel("orbits")
    ax_max.set_ylabel(r"$\alpha$")
    ax_max.set_title(r"$\alpha_M$")
    ax_max.set_ylim([-0.1,1])
    ax_max.set_xlim([40, 60])

    ax_freyn = fig.add_subplot(gs[1, 0])
    ax_freyn.set_xlabel("orbits")
    if B_field_focus:
        ax_freyn.set_ylabel(r"P_B")
        ax_freyn.set_title("Magnetic Pressure")
    else:
        ax_freyn.set_ylabel(r"$\alpha$")
        ax_freyn.set_title(r"$\rho v_r v_{\phi}$")

    ax_treyn = fig.add_subplot(gs[0, 1])
    ax_treyn.set_xlabel("orbits")
    if B_field_focus:
        ax_treyn.set_ylabel(r"$B_z$")
        ax_treyn.set_title(r"$B_z$")
    else:
        ax_treyn.set_ylabel(r"$\alpha$")
        ax_treyn.set_title(r"$\rho \left(v_r - \left<v_r\right>\right) \left(v_{\phi} - \left<v_{\phi}\right>\right)$")

    if not B_field_focus:
        ax_creyn = fig.add_subplot(gs[1, 1])
        ax_creyn.set_xlabel("orbits")
        ax_creyn.set_ylabel(r"$\alpha$")
        ax_creyn.set_title(r"$\rho v_r \left(v_{\phi} - \left<v_{\phi}\right>\right)$")

    for d in range(len(dnames)):
        if final_points[d] == 0:
            ax_max.plot(orbits[:fnum], alphas[d, :fnum, 0], f"C{d}-", label=dnames[d])
            ax_treyn.plot(orbits[:fnum], alphas[d, :fnum, 1], f"C{d}-", label=dnames[d])
            ax_freyn.plot(orbits[:fnum], alphas[d, :fnum, 2], f"C{d}-", label=dnames[d])
            if not B_field_focus:
                ax_creyn.plot(orbits[:fnum], alphas[d, :fnum, 3], f"C{d}-", label=dnames[d])
        else:
            ax_max.plot(orbits[:int(final_points[d])], alphas[d, :int(final_points[d]), 0], f"C{d}-", label=dnames[d])
            ax_treyn.plot(orbits[:int(final_points[d])], alphas[d, :int(final_points[d]), 1], f"C{d}-", label=dnames[d])
            ax_freyn.plot(orbits[:int(final_points[d])], alphas[d, :int(final_points[d]), 2], f"C{d}-", label=dnames[d])
            if not B_field_focus:
                ax_creyn.plot(orbits[:int(final_points[d])], alphas[d, :int(final_points[d]), 3], f"C{d}-", label=dnames[d])
    
    plt.legend()
    plt.subplots_adjust(top=(1-0.01*(16/vert)))
    if B_field_focus:
        title = "Mass Weighted Average quantities"
    else:
        title = r"Mass Weighted Average $\alpha$s"
    if res:
        title += " at Resonant Radius"
    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig("%s/%s.png" % (savedir, aname))
    plt.close()

