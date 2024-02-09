from .roots.athena_analysis import *
from branches.roots.misc_func import *
import pickle
import scipy.fft

logging.basicConfig(filename=file.logs_loc+"/waves.log", encoding='utf-8', level=logging.INFO)

def fourier_waves(dname, fnum, m_range):
    aname = "_waves" #a for analysis
    data_location = file.data_loc + dname
    grid_type=file.grid_types[dname]
    savedir = file.savedir + dname + "/" + dname + aname
    mkdir_if_not_exist(savedir)
    filename = "%s/disk.out1.%05d.athdf" % (data_location, fnum)
    aa = Athena_Analysis(filename=filename, grid_type=grid_type)
    m_range = np.arange(m_range[0], m_range[1]+1, 1)

    aa.native_grid()
    orbit = (aa.time / sim.binary_period)
    aa.get_primaries(get_rho=True)

    vert = 1
    horz = 1
    gs = gridspec.GridSpec(vert, horz)
    fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
    ax = fig.add_subplot(gs[0, 0])

    maxs = np.zeros(len(m_range))
    mins = np.zeros(len(m_range))

    for m in m_range:
        disk_rho = aa.integrate(aa.rho, "vert")
        fourier_wave_base = np.zeros(disk_rho.shape)
        differential = np.zeros(disk_rho.shape)
        for n in range(disk_rho.shape[1]):
            fourier_wave_base[:,n] = np.arange(0, 2*np.pi, (2*np.pi/disk_rho.shape[0]))
            differential[:,n] = aa.possible_dphi_primitive
        fourier_wavec = np.cos(m*fourier_wave_base)
        fourier_waves = np.sin(m*fourier_wave_base)
        fourier_coefficients = np.sqrt((np.sum((fourier_wavec*disk_rho*differential), axis=0))**2
             + (np.sum((fourier_waves*disk_rho*differential), axis=0))**2) / (2*np.pi)

        ax.plot(aa.possible_r, fourier_coefficients, f"C{m}-", label=f"m = {m}")
        maxs[m-m_range[0]] = max(fourier_coefficients)
        mins[m-m_range[0]] = min(fourier_coefficients)

    true_max = max(maxs)
    true_min = min(mins)

    ax.set_ylim([-15,35])
    true_max = 35
    true_min = -15

    vertical_distribution = np.arange(true_min, true_max, 1)

    ax.plot(np.full(len(vertical_distribution), sim.three_one_res), vertical_distribution, "C3--")

    ax.set_xlabel("r")
    ax.set_ylabel("fourier coefficient")
    ax.set_title(f"orbit: {orbit:.2f}")
    ax.legend()
    plt.savefig("%s/%s%s%05d.png" % (savedir, dname, aname, fnum))
    plt.close()

def fourier_waves_loop(dname, fnum_range, file_spacing, m_range, grid_type):
    now = datetime.now()
    fnum_range = np.arange(fnum_range[0], fnum_range[-1]+1, file_spacing)
    for fnum in fnum_range:
        logging.info(datetime.now()-now)
        now = datetime.now()
        logging.info("fnum = %d" % fnum)
        fourier_waves(dname, fnum, m_range, grid_type)

class Fourier_Waves:
    def __init__(self, dname, file_spacing=None):
        self.dname = dname
        self.aname = "_waves"
        if file_spacing is None:
            self.file_spacing = find_file_spacing(dname)
        else:
            self.file_spacing = file_spacing
        self.data_location = file.data_loc + self.dname
        self.grid_type=file.grid_types[self.dname]
        self.savedir = file.savedir + self.dname + "/" + self.dname + self.aname
        mkdir_if_not_exist(self.savedir)
        self.pickldir = self.savedir + "/pickles"
        mkdir_if_not_exist(self.pickldir)

    def collect_density(self, start_fnum=0, start_orbit=None, orbits=3):
        if start_orbit is not None:
            start_fnum = start_orbit * sim.filenums_per_orbit
        aname = "density"
        
        #finding file range
        now = datetime.now()
        max_fnum = self.file_spacing * math.ceil((orbits*sim.filenums_per_orbit+start_fnum)/self.file_spacing)
        if not os.path.exists("%s/disk.out1.%05d.athdf" % (self.data_location, max_fnum)):
            raise("Not enough data left, start further back")
        fnums = np.arange(start_fnum, max_fnum+self.file_spacing, self.file_spacing)
        file_count = ((max_fnum-start_fnum) / self.file_spacing) + 1
        logging.info(f"Found data, range: [{start_fnum}, {max_fnum}], file_count: {file_count}") #purposefully didn't make file count an int yet, so you can verify it is an int in log
        file_count = int(file_count)

        logging.info("Beginning collection...")
        for i, fnum in enumerate(fnums):
            logging.info(datetime.now()-now)
            now = datetime.now()
            logging.info("fnum = %d" % fnum)
            
            filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
            aa = Athena_Analysis(filename=filename, grid_type=self.grid_type) 

            if i == 0:
                aa._axes()
                data = np.zeros((file_count, len(aa.possible_phi), len(aa.possible_r)))
                normalization = aa.integrate(1, "z")
            
            if i != ((fnum-start_fnum) / self.file_spacing):
                raise(f"Something is very wrong, fnum / filespacing = {(fnum-start_fnum) / self.file_spacing}]")
            aa.get_primaries(get_rho=True)
            data[i] = aa.integrate(aa.rho, 'z') / normalization

            if (fnum == max_fnum):
                logging.info("pickling up to %s" % fnum)
                with open("%s/%s%s_pickle_%05d-%05d.dat" % (self.pickldir, self.dname, aname, start_fnum, fnum), "wb") as pickle_file:
                    pickle.dump({"data": data, "fnums": fnums, "phis": aa.possible_phi, "rs": aa.possible_r,
                        "dphi": aa.possible_dphi_primitive, "dr": aa.possible_dr_primitive}, pickle_file)
                self.pickle_loc = "%s/%s_pickle_%05d-%05d.dat" % (self.pickldir, self.dname, start_fnum, fnum)

    def collect_momentum_density(self, start_fnum=0, start_orbit=None, orbits=3):
        if start_orbit is not None:
            start_fnum = start_orbit * sim.filenums_per_orbit
        aname = "momentum_density"
        
        #finding file range
        now = datetime.now()
        max_fnum = self.file_spacing * math.ceil((orbits*sim.filenums_per_orbit+start_fnum)/self.file_spacing)
        if not os.path.exists("%s/disk.out1.%05d.athdf" % (self.data_location, max_fnum)):
            raise("Not enough data left, start further back")
        fnums = np.arange(start_fnum, max_fnum+self.file_spacing, self.file_spacing)
        file_count = ((max_fnum-start_fnum) / self.file_spacing) + 1
        logging.info(f"Found data, range: [{start_fnum}, {max_fnum}], file_count: {file_count}") #purposefully didn't make file count an int yet, so you can verify it is an int in log
        file_count = int(file_count)

        logging.info("Beginning collection...")
        for i, fnum in enumerate(fnums):
            logging.info(datetime.now()-now)
            now = datetime.now()
            logging.info("fnum = %d" % fnum)
            
            filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
            aa = Athena_Analysis(filename=filename, grid_type=self.grid_type) 

            if i == 0:
                aa._axes()
                data = np.zeros((3, file_count, len(aa.possible_phi), len(aa.possible_r)))
                normalization = aa.integrate(1, "z")
            
            if i != ((fnum-start_fnum) / self.file_spacing):
                raise(f"Something is very wrong, fnum / filespacing = {(fnum-start_fnum) / self.file_spacing}]")
            aa.get_primaries(get_rho=True, get_vel_r=True, get_vel_phi=True, get_vel_z=True)
            data[0, i] = aa.integrate(aa.rho * aa.vel_z, 'z') / normalization
            data[1, i] = aa.integrate(aa.rho * aa.vel_phi, 'z') / normalization
            data[2, i] = aa.integrate(aa.rho * aa.vel_r, 'z') / normalization

            if (fnum == max_fnum):
                logging.info("pickling up to %s" % fnum)
                with open("%s/%s%s_pickle_%05d-%05d.dat" % (self.pickldir, self.dname, aname, start_fnum, fnum), "wb") as pickle_file:
                    pickle.dump({"data": data, "fnums": fnums, "phis": aa.possible_phi, "rs": aa.possible_r,
                        "dphi": aa.possible_dphi_primitive, "dr": aa.possible_dr_primitive}, pickle_file)
                self.pickle_loc = "%s/%s_pickle_%05d-%05d.dat" % (self.pickldir, self.dname, start_fnum, fnum)

    def collect_rholrl(self, start_fnum=0, start_orbit=None, orbits=3):
            if start_orbit is not None:
                start_fnum = start_orbit * sim.filenums_per_orbit
            aname = "rholrl"
            
            #finding file range
            now = datetime.now()
            max_fnum = self.file_spacing * math.ceil((orbits*sim.filenums_per_orbit+start_fnum)/self.file_spacing)
            if not os.path.exists("%s/disk.out1.%05d.athdf" % (self.data_location, max_fnum)):
                raise("Not enough data left, start further back")
            fnums = np.arange(start_fnum, max_fnum+self.file_spacing, self.file_spacing)
            file_count = ((max_fnum-start_fnum) / self.file_spacing) + 1
            logging.info(f"Found data, range: [{start_fnum}, {max_fnum}], file_count: {file_count}") #purposefully didn't make file count an int yet, so you can verify it is an int in log
            file_count = int(file_count)

            logging.info("Beginning collection...")
            for i, fnum in enumerate(fnums):
                logging.info(datetime.now()-now)
                now = datetime.now()
                logging.info("fnum = %d" % fnum)
                
                filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
                aa = Athena_Analysis(filename=filename, grid_type=self.grid_type) 

                if i == 0:
                    aa._axes()
                    data = np.zeros((3, file_count, len(aa.possible_phi), len(aa.possible_r)))
                    normalization = aa.integrate(1, "z")
                
                if i != ((fnum-start_fnum) / self.file_spacing):
                    raise(f"Something is very wrong, fnum / filespacing = {(fnum-start_fnum) / self.file_spacing}]")
                aa.get_primaries(get_rho=True, get_vel_r=True, get_vel_phi=True, get_vel_z=True)
                
                lrl_native, lrl_cart = aa.get_lrl(components = True)
                rhoxlrl = np.array([aa.rho*lrl_cart[0], aa.rho*lrl_cart[1], np.zeros(aa.array_size)])
                data[0, i] = aa.integrate(rhoxlrl[0], "z") / normalization
                data[1, i] = aa.integrate(rhoxlrl[1], "z") / normalization
                data[2, i] = aa.integrate(rhoxlrl[2], "z") / normalization

                if (fnum == max_fnum):
                    logging.info("pickling up to %s" % fnum)
                    with open("%s/%s%s_pickle_%05d-%05d.dat" % (self.pickldir, self.dname, aname, start_fnum, fnum), "wb") as pickle_file:
                        pickle.dump({"data": data, "fnums": fnums, "phis": aa.possible_phi, "rs": aa.possible_r,
                            "dphi": aa.possible_dphi_primitive, "dr": aa.possible_dr_primitive}, pickle_file)
                    self.pickle_loc = "%s/%s_pickle_%05d-%05d.dat" % (self.pickldir, self.dname, start_fnum, fnum)
    
    def plot_density(self, fnum_range=None, r_slicepoint=10, bounds=[0, 20], log=False):
        aname = "density"
        if fnum_range is None:
            pickl_loc = self.pickle_loc
        else:
            start_fnum = fnum_range[0]
            max_fnum = fnum_range[1]
            pickl_loc = "%s/%s%s_pickle_%05d-%05d.dat" % (self.pickldir, self.dname, aname, start_fnum, max_fnum)

        with open(pickl_loc, "rb") as pickle_file:
            pickl = pickle.load(pickle_file)

        phis = pickl["phis"]
        dphi = pickl["dphi"]
        phi_ax = pickl["phis"] / (np.pi)
        times = pickl["fnums"]
        dt = 1 / sim.filenums_per_orbit
        time_ax = pickl["fnums"] / sim.filenums_per_orbit
        r_ax = pickl["rs"]
        r_idx = np.argmin(abs(r_ax - r_slicepoint))
        data = pickl["data"]
        
        phi_fourier = np.zeros((len(time_ax), len(r_ax)), dtype = 'complex_')
        for t in range(len(times)):
            for r in range(len(r_ax)):
                phi_fourier[t, r] = np.sum(data[t, :, r] * np.exp(2j*phis) * dphi) / 2*np.pi

        t_fourier = np.zeros((len(phi_ax), len(r_ax)), dtype = 'complex_')
        for p in range(len(phis)):
            for r in range(len(r_ax)):
                t_fourier[p, r] = np.sum(data[:, p, r] * np.exp(-1j*(2*np.pi)*time_ax) * dt) / (time_ax[-1]-time_ax[0])

        phi_t_fourier = np.zeros(len(r_ax), dtype = 'complex_')
        for r in range(len(r_ax)):
            phi_t_fourier[r] = np.sum(phi_fourier[:, r] * np.exp(-1j*(2*np.pi)*time_ax) * dt) / (time_ax[-1]-time_ax[0])
        
        vert = 2
        horz = 4
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
        ax = fig.add_subplot(gs[0, 0])
        ax_p = fig.add_subplot(gs[1, 0])
        ax_t = fig.add_subplot(gs[1, 1])
        ax_pt = fig.add_subplot(gs[0, 2])
        ax_pta = fig.add_subplot(gs[1, 2])
        ax_re = fig.add_subplot(gs[0, 3])
        ax_im = fig.add_subplot(gs[1, 3])


        clip=True
        if log == True:
            norm = colors.LogNorm(bounds[0], bounds[1], clip)
            normt = colors.LogNorm(bounds[0], bounds[1]/2, clip)
            normrho = colors.LogNorm(bounds[0], bounds[1]/2, clip)
        else:
            norm = colors.Normalize(bounds[0], bounds[1], clip)
            normt = colors.Normalize(bounds[0], bounds[1]/2, clip)
            normrho = colors.Normalize(bounds[0], bounds[1]/2, clip)

        im = ax.pcolormesh(phi_ax, time_ax, data[:,:,r_idx], cmap="viridis", norm=normrho, clip_on=True, shading="auto")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        ax.grid(False)
        cbar = plt.colorbar(im, cax=cax)
        cbar.ax.tick_params()
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
        ax.set_xlabel(r"$\phi$")
        ax.set_ylabel("Binary Orbit")
        ax.set_title(r"$\rho$"+f" r={r_slicepoint:.3}")

        im_p = ax_p.pcolormesh(r_ax, time_ax, np.abs(phi_fourier), cmap="viridis", norm=norm, clip_on=True, shading="auto")
        divider = make_axes_locatable(ax_p)
        cax_p = divider.append_axes("right", size="5%", pad=0.05)
        ax_p.grid(False)
        cbar = plt.colorbar(im_p, cax=cax_p)
        cbar.ax.tick_params()
        ax_p.set_xlabel(r"r")
        ax_p.set_ylabel("Binary Orbit")
        ax_p.set_title(r"fourier coefficient magnitude $e^{2i\phi}$")

        im_t = ax_t.pcolormesh(r_ax, phi_ax, np.abs(t_fourier), cmap="viridis", norm=normt, clip_on=True, shading="auto")
        divider = make_axes_locatable(ax_t)
        cax_t = divider.append_axes("right", size="5%", pad=0.05)
        ax_t.grid(False)
        cbar = plt.colorbar(im_t, cax=cax_t)
        cbar.ax.tick_params()
        ax_t.set_xlabel(r"r")
        ax_t.set_ylabel(r"$phi$")
        ax_t.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
        ax_t.set_title(r"fourier coefficient magnitude $e^{-i\Omega_b t}$")

        ax_pt.plot(r_ax, np.abs(phi_t_fourier))
        ax_pt.axvline(x= sim.three_one_res, color = "red", linestyle="--")
        ax_pt.set_xlabel("r")
        ax_pt.set_title(r"Magnitude of Fourier: $e^{i\left(2\phi - \Omega_b t\right)}$")

        ax_pta.plot(r_ax, np.angle(phi_t_fourier))
        ax_pta.axvline(x= sim.three_one_res, color = "red", linestyle="--")
        ax_pta.set_xlabel("r")
        ax_pta.set_title(r"Argument of Fourier: $e^{i\left(2\phi - \Omega_b t\right)}$")

        ax_re.plot(r_ax, np.real(phi_t_fourier))
        ax_re.axvline(x= sim.three_one_res, color = "red", linestyle="--")
        ax_re.set_xlabel("r")
        ax_re.set_title(r"Real Part of Fourier: $e^{i\left(2\phi - \Omega_b t\right)}$")

        ax_im.plot(r_ax, np.imag(phi_t_fourier))
        ax_im.axvline(x= sim.three_one_res, color = "red", linestyle="--")
        ax_im.set_xlabel("r")
        ax_im.set_title(r"Imaginary Part of Fourier: $e^{i\left(2\phi - \Omega_b t\right)}$")

        plt.tight_layout()
        plt.savefig("%s/%s%s%05d-%05d.png" % (self.savedir, self.dname, self.aname, times[0], times[-1]))
        plt.close()

    def plot_momentum_density(self, fnum_range=None, r_slicepoint=10, bounds=[0, 20], log=False):
        aname = "momentum_density"
        if fnum_range is None:
            pickl_loc = self.pickle_loc
        else:
            start_fnum = fnum_range[0]
            max_fnum = fnum_range[1]
            pickl_loc = "%s/%s%s_pickle_%05d-%05d.dat" % (self.pickldir, self.dname, aname, start_fnum, max_fnum)

        with open(pickl_loc, "rb") as pickle_file:
            pickl = pickle.load(pickle_file)

        phis = pickl["phis"]
        dphi = pickl["dphi"]
        phi_ax = pickl["phis"] / (np.pi)
        times = pickl["fnums"]
        dt = 1 / sim.filenums_per_orbit
        time_ax = pickl["fnums"] / sim.filenums_per_orbit
        r_ax = pickl["rs"]
        r_idx = np.argmin(abs(r_ax - r_slicepoint))
        data = pickl["data"]
        
        coords=["z", "phi", "r"]
        for idx in range(3):
            sname = coords[idx]

            phi_fourier = np.zeros((len(time_ax), len(r_ax)), dtype = 'complex_')
            for t in range(len(times)):
                for r in range(len(r_ax)):
                    phi_fourier[t, r] = np.sum(data[idx, t, :, r] * np.exp(2j*phis) * dphi) / 2*np.pi

            t_fourier = np.zeros((len(phi_ax), len(r_ax)), dtype = 'complex_')
            for p in range(len(phis)):
                for r in range(len(r_ax)):
                    t_fourier[p, r] = np.sum(data[idx, :, p, r] * np.exp(-1j*(2*np.pi)*time_ax) * dt) / (time_ax[-1]-time_ax[0])

            phi_t_fourier = np.zeros(len(r_ax), dtype = 'complex_')
            for r in range(len(r_ax)):
                phi_t_fourier[r] = np.sum(phi_fourier[:, r] * np.exp(-1j*(2*np.pi)*time_ax) * dt) / (time_ax[-1]-time_ax[0])
            
            vert = 2
            horz = 4
            gs = gridspec.GridSpec(vert, horz)
            fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
            ax = fig.add_subplot(gs[0, 0])
            ax_p = fig.add_subplot(gs[1, 0])
            ax_t = fig.add_subplot(gs[1, 1])
            ax_pt = fig.add_subplot(gs[0, 2])
            ax_pta = fig.add_subplot(gs[1, 2])
            ax_re = fig.add_subplot(gs[0, 3])
            ax_im = fig.add_subplot(gs[1, 3])


            clip=True
            if log == True:
                norm = colors.LogNorm(bounds[0], bounds[1], clip)
                normt = colors.LogNorm(bounds[0], bounds[1]/2, clip)
                normrho = colors.LogNorm(bounds[0], bounds[1]/2, clip)
            else:
                norm = colors.Normalize(bounds[0], bounds[1], clip)
                normt = colors.Normalize(bounds[0], bounds[1]/2, clip)
                normrho = colors.Normalize(bounds[0], bounds[1]/2, clip)

            im = ax.pcolormesh(phi_ax, time_ax, data[idx,:,:,r_idx], cmap="viridis", norm=normrho, clip_on=True, shading="auto")
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            ax.grid(False)
            cbar = plt.colorbar(im, cax=cax)
            cbar.ax.tick_params()
            ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
            ax.set_xlabel(r"$\phi$")
            ax.set_ylabel("Binary Orbit")
            ax.set_title(f"sliced at r={r_slicepoint:.3}")

            im_p = ax_p.pcolormesh(r_ax, time_ax, np.abs(phi_fourier), cmap="viridis", norm=norm, clip_on=True, shading="auto")
            divider = make_axes_locatable(ax_p)
            cax_p = divider.append_axes("right", size="5%", pad=0.05)
            ax_p.grid(False)
            cbar = plt.colorbar(im_p, cax=cax_p)
            cbar.ax.tick_params()
            ax_p.set_xlabel(r"r")
            ax_p.set_ylabel("Binary Orbit")
            ax_p.set_title(r"fourier coefficient magnitude $e^{2i\phi}$")

            im_t = ax_t.pcolormesh(r_ax, phi_ax, np.abs(t_fourier), cmap="viridis", norm=normt, clip_on=True, shading="auto")
            divider = make_axes_locatable(ax_t)
            cax_t = divider.append_axes("right", size="5%", pad=0.05)
            ax_t.grid(False)
            cbar = plt.colorbar(im_t, cax=cax_t)
            cbar.ax.tick_params()
            ax_t.set_xlabel(r"r")
            ax_t.set_ylabel(r"$phi$")
            ax_t.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
            ax_t.set_title(r"fourier coefficient magnitude $e^{-i\Omega_b t}$")

            ax_pt.plot(r_ax, np.abs(phi_t_fourier))
            ax_pt.axvline(x= sim.three_one_res, color = "red", linestyle="--")
            ax_pt.set_xlabel("r")
            ax_pt.set_title(r"Magnitude of Fourier: $e^{i\left(2\phi - \Omega_b t\right)}$")

            ax_pta.plot(r_ax, np.angle(phi_t_fourier))
            ax_pta.axvline(x= sim.three_one_res, color = "red", linestyle="--")
            ax_pta.set_xlabel("r")
            ax_pta.set_title(r"Argument of Fourier: $e^{i\left(2\phi - \Omega_b t\right)}$")

            ax_re.plot(r_ax, np.real(phi_t_fourier))
            ax_re.axvline(x= sim.three_one_res, color = "red", linestyle="--")
            ax_re.set_xlabel("r")
            ax_re.set_title(r"Real Part of Fourier: $e^{i\left(2\phi - \Omega_b t\right)}$")

            ax_im.plot(r_ax, np.imag(phi_t_fourier))
            ax_im.axvline(x= sim.three_one_res, color = "red", linestyle="--")
            ax_im.set_xlabel("r")
            ax_im.set_title(r"Imaginary Part of Fourier: $e^{i\left(2\phi - \Omega_b t\right)}$")

            plt.subplots_adjust(top=(1-0.01*(16/vert)))
            fig.suptitle(f"{coords[idx]} momentum density")
            plt.tight_layout()
            plt.savefig("%s/%s%s%05d-%05d%s.png" % (self.savedir, self.dname, self.aname, times[0], times[-1], sname))
            plt.close()

    def tidal_waves(self, fnum_range=None):
        aname = "tidal_waves"
        if fnum_range is None:
            pickl_loc = self.pickle_loc
        else:
            start_fnum = fnum_range[0]
            max_fnum = fnum_range[1]
            pickl_loc = "%s/%s%s_pickle_%05d-%05d.dat" % (self.pickldir, self.dname, "momentum_density", start_fnum, max_fnum)

        with open(pickl_loc, "rb") as pickle_file:
            pickl = pickle.load(pickle_file)

        phis = pickl["phis"]
        dphi = pickl["dphi"]
        phi_ax = pickl["phis"] / (np.pi)
        times = pickl["fnums"]
        dt = 1 / sim.filenums_per_orbit
        time_ax = pickl["fnums"] / sim.filenums_per_orbit
        r_ax = pickl["rs"]
        data = pickl["data"]

        filename = "%s/disk.out1.%05d.athdf" % (self.data_location, int((fnum_range[1]-fnum_range[0])/2))
        ref = Athena_Analysis(filename=filename, grid_type=self.grid_type)
        
        ref.get_potentials(get_companion_grav=True, get_accel=True)
        ref.get_primaries(get_rho=True)
        M_fluid = ref.integrate(ref.rho, "all")
        tidal_field = ref.accel_pot + ref.companion_grav_pot

        tides = ref.integrate(tidal_field, "z") / ref.integrate(1, "z")

        wave_modes = [2]
        C_m = np.zeros((len(wave_modes), len(r_ax)))
        for i, m in enumerate(wave_modes):
            
            tides_phi_fourier = np.zeros((len(time_ax), len(r_ax)), dtype = 'complex_')
            pr_phi_fourier = np.zeros((len(time_ax), len(r_ax)), dtype = 'complex_')
            pphi_phi_fourier = np.zeros((len(time_ax), len(r_ax)), dtype = 'complex_')
            for t in range(len(times)):
                for r in range(len(r_ax)):
                    current_phis = (phis + sim.orbital_Omega * times[t])
                    tides_phi_fourier[t, r] = np.sum(tides[:, r] * np.exp(1j*(m)*current_phis * dphi)) / 2*np.pi
                    pr_phi_fourier[t, r] = np.sum(data[2, t, :, r] * np.exp(-1j*(m+1)*current_phis) * dphi) / 2*np.pi
                    pphi_phi_fourier[t, r] = np.sum(data[1, t, :, r] * np.exp(-1j*(m+1)*current_phis) * dphi) / 2*np.pi
        
            tid_m = np.zeros(len(r_ax), dtype = 'complex_')
            pr_m = np.zeros(len(r_ax), dtype = 'complex_')
            pphi_m = np.zeros(len(r_ax), dtype = 'complex_')
            for r in range(len(r_ax)):
                tid_m[r] = np.sum(tides_phi_fourier[:, r] * np.exp(-1j*(m)*sim.orbital_Omega*time_ax) * dt) / (time_ax[-1]-time_ax[0])
                pr_m[r] = np.sum(pr_phi_fourier[:, r] * np.exp(-1j*(m)*sim.orbital_Omega*time_ax) * dt) / (time_ax[-1]-time_ax[0])
                pphi_m[r] = np.sum(pphi_phi_fourier[:, r] * np.exp(-1j*(m)*sim.orbital_Omega*time_ax) * dt) / (time_ax[-1]-time_ax[0])

            d_tid_m_dr = np.zeros(len(tid_m))
            d_tid_m_dr[1:] = (tid_m[1:] - tid_m[:-1]) / (r_ax[1:] - r_ax[:-1])
            d_tid_m_dr[0] = d_tid_m_dr[1] - ((d_tid_m_dr[2]-d_tid_m_dr[1]) * ((r_ax[1]-r_ax[0])/(r_ax[2]-r_ax[1])))

            prefactor = (2*np.pi)*(2*np.pi)/(2*sim.gm1*M_fluid*sim.orbital_Omega)
            expression = -1*m*tid_m*(2j*pphi_m + pr_m) + 1j*r_ax*pphi_m*d_tid_m_dr
            C_m[i] = np.real_if_close(prefactor*(expression + np.conjugate(expression)))

            
        vert = 1
        horz = 1
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
        ax = fig.add_subplot(gs[0, 0])
        for i, C in enumerate(C_m):
            ax.plot(r_ax, C, f"C{i}-", label=f"({m[i]},{m[i]+1})")
        ax.axvline(x= sim.three_one_res, color = "red", linestyle="--")
        ax.set_xlabel("r")
        ax.set_ylabel("Eccent Contribution")
        ax.legend(loc="upper left")
        
        plt.tight_layout()
        plt.savefig("%s/%s%s%05d-%05d.png" % (self.savedir, self.dname, self.aname, times[0], times[-1]))
        plt.close()