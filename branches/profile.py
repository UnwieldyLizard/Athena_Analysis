#Eventually needs to contain rewrite of mprofile code
from .roots.athena_analysis import *
from .roots.misc_func import *
import pickle


logging.basicConfig(filename=file.logs_loc+"/profile.log", encoding='utf-8', level=logging.INFO)

class Profile():
    def __init__(self, dname):
        self.dname = dname
        self.aname = "_profile" #a for analysis
        self.data_location = file.data_loc + dname
        self.grid_type = file.grid_types[dname]
        self.is_MHD = file.MHD[dname]
        self.savedir = file.savedir + dname + "/" + dname + self.aname
        mkdir_if_not_exist(self.savedir)

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

        aa.midplane_colorplot(plasma_beta, ax, vbound=[1e-3,1e5], log=True)
        aa.midplane_colorplot(plasma_beta, ax_vert, vbound=[1e-3,1e5], log=True, slicetype='y')

        orbit = (aa.time / sim.binary_period)
        plt.subplots_adjust(top=(1-0.01*(16/vert)))
        fig.suptitle(f"orbit: {orbit:.2f}")
        plt.savefig("%s/%s%s%05d.png" % (self.savedir, self.dname, aname, fnum))
        plt.close()

def compare_beta(dnames, fnum_range, inverse=False):
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

    orbits = np.zeros([(fnum_range[1] - fnum_range[0] + offset)])
    betas = np.zeros([len(dnames), (fnum_range[1] - fnum_range[0] + offset)])

    if dnames != data["dnames"]:
        raise("this functionality doesn't exist yet, please use consistent dnames")
    if load_point != 0:
        orbits[:offset] = data["orbits"]
        betas[:, :offset] = data["betas"][:,:offset]

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
            aa = Athena_Analysis(filename=filename, grid_type=file.grid_types[dnames[d]])
            
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

            ax.plot(orbits[:i+1], betas[d, :i+1], f"C{d}-", label=dnames[d])

            with open("%s/pickles/pickle_%05d.dat" % (savedir, f), "wb") as pickle_file:
                pickle.dump({
                    "inverse": inverse,
                    "dnames": dnames,
                    "orbits": orbits[:+i],
                    "betas": betas[:,:+1],
                }, pickle_file)

        plt.legend()
        plt.savefig("%s/%s%05d.png" % (savedir, aname, f))
        plt.close()

def compare_alpha(dnames, fnum_range):
    aname = "alpha"
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

    orbits = np.zeros([(fnum_range[1] - fnum_range[0] + offset)])
    alphas = np.zeros([len(dnames), (fnum_range[1] - fnum_range[0] + offset)])
    
    if load_point != 0:
        if dnames != data["dnames"]:
            raise("this functionality doesn't exist yet, please use consistent dnames")
        orbits[:offset] = data["orbits"]
        alphas[:, :offset] = data["alphas"][:,:offset]

    for j, f in enumerate(np.arange(fnum_range[0], fnum_range[1], file_spacing)):
        i = j + offset
        vert = 1
        horz = 1
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
        ax = fig.add_subplot(gs[0, 0])
        ax.set_xlabel("orbits")
        #ax.set_yscale("log")
        ax.set_ylabel(r"$\alpha$")
        ax.set_title(r"Mass Weighted Average $\alpha$")
        for d in range(len(dnames)):
            filename = "%s/disk.out1.%05d.athdf" % (file.data_loc + dnames[d], f)
            aa = Athena_Analysis(filename=filename, grid_type=file.grid_types[dnames[d]])
            
            aa.get_primaries(get_press=True, get_rho=True)
            aa.get_Bfields()
            B_press = (aa.B_z**2 + aa.B_phi**2 + aa.B_r**2)/2

            averaged_Maxwell = aa.integrate(aa.rho*(-2/3)*(aa.B_r*aa.B_phi)/(aa.press+B_press), "All") / aa.integrate(aa.rho, "All")
            
            alphas[d, i] = averaged_Maxwell
            orbits[i] = aa.time / sim.binary_period

            ax.plot(orbits[:i+1], alphas[d, :i+1], f"C{d}-", label=dnames[d])

            print(alphas[d,:i+1])

            with open("%s/pickles/pickle_%05d.dat" % (savedir, f), "wb") as pickle_file:
                pickle.dump({
                    "dnames": dnames,
                    "orbits": orbits[:i+1],
                    "alphas": alphas[:,:i+1],
                }, pickle_file)

        plt.legend()
        plt.savefig("%s/%s%05d.png" % (savedir, aname, f))
        plt.close()
    