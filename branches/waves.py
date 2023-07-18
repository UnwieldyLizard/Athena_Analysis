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
    def __init__(self, dname=None):
        self.dname = dname
        self.aname = "_waves"
    def collect_data(self, dname=None, start_fnum=0, orbits=3):
        if dname is None:
            dname = self.dname
        data_location = file.data_loc + dname
        grid_type=file.grid_types[dname]
        savedir = file.savedir + dname + "/" + dname + self.aname
        mkdir_if_not_exist(savedir)
        pickldir = savedir + "/pickles"
        mkdir_if_not_exist(pickldir)

        #finding file range
        now = datetime.now()
        logging.info("Finding files...")
        max_fnum = 1
        file_spacing = None
        while file_spacing is None:
            while not os.path.exists("%s/disk.out1.%05d.athdf" % (data_location, max_fnum)):
                max_fnum += 1
            file_spacing = max_fnum
        logging.info(f"Found file spacing: {file_spacing}")
        max_fnum = file_spacing * math.ceil((orbits*sim.filenums_per_orbit+start_fnum)/file_spacing)
        if not os.path.exists("%s/disk.out1.%05d.athdf" % (data_location, max_fnum)):
            raise("Not enough data left, start further back")
        fnums = np.arange(start_fnum, max_fnum+file_spacing, file_spacing)
        file_count = ((max_fnum-start_fnum) / file_spacing) + 1
        logging.info(f"Found data, range: [{start_fnum}, {max_fnum}], file_count: {file_count}") #purposefully didn't make file count an int yet, so you can verify it is an int in log
        file_count = int(file_count)

        logging.info("Beginning collection...")
        for i, fnum in enumerate(fnums):
            logging.info(datetime.now()-now)
            now = datetime.now()
            logging.info("fnum = %d" % fnum)
            
            filename = "%s/disk.out1.%05d.athdf" % (data_location, fnum)
            aa = Athena_Analysis(filename=filename, grid_type=grid_type) 

            if i == 0:
                aa._axes()
                data = np.zeros((file_count, len(aa.possible_phi), len(aa.possible_r)))
                normalization = aa.integrate(1, "z")
            
            if i != ((fnum-start_fnum) / file_spacing):
                raise(f"Something is very wrong, fnum / filespacing = {(fnum-start_fnum) / file_spacing}]")
            aa.get_primaries(get_rho=True)
            data[i] = aa.integrate(aa.rho, 'z') / normalization

            if (fnum == max_fnum):
                logging.info("pickling up to %s" % fnum)
                with open("%s/%s_pickle_%05d-%05d.dat" % (pickldir, dname, start_fnum, fnum), "wb") as pickle_file:
                    pickle.dump({"data": data, "fnums": fnums, "phis": aa.possible_phi, "rs": aa.possible_r,
                        "dphi": aa.possible_dphi_primitive, "dr": aa.possible_dr_primitive}, pickle_file)
                #remembering pickle loc
                self.pickle_loc = "%s/%s_pickle_%05d-%05d.dat" % (pickldir, dname, start_fnum, fnum)
    
    def plot(self, dname=None, fnum_range=None, r_slicepoint=10, bounds=[0, 20], log=False):
        if dname is None:
            dname = self.dname
        data_location = file.data_loc + dname
        grid_type=file.grid_types[dname]
        savedir = file.savedir + dname + "/" + dname + self.aname
        mkdir_if_not_exist(savedir)
        pickldir = savedir + "/pickles"
        mkdir_if_not_exist(pickldir)
        if fnum_range is None:
            pickl_loc = self.pickle_loc
        else:
            start_fnum = fnum_range[0]
            max_fnum = fnum_range[1]
            pickl_loc = "%s/%s_pickle_%05d-%05d.dat" % (pickldir, dname, start_fnum, max_fnum)

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
        horz = 2
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
        ax = fig.add_subplot(gs[0, 0])
        ax_p = fig.add_subplot(gs[1, 0])
        ax_t = fig.add_subplot(gs[1, 1])
        ax_pt = fig.add_subplot(gs[0, 1])

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
        ax.set_title(r"$\rho$"+f" r={r_slicepoint}")

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

        ax_pt.plot(r_ax, abs(phi_t_fourier))
        ax_pt.set_xlabel("r")
        ax_pt.set_title(r"fourier coefficient magnitude $e^{i\left(2\phi - \Omega_b t\right)}$")

        plt.tight_layout()
        plt.savefig("%s/%s%s%05d-%05d.png" % (savedir, dname, self.aname, times[0], times[-1]))
        plt.close()