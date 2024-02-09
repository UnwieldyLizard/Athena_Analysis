from .roots.athena_analysis import *
from .roots.misc_func import *
import pickle

logging.basicConfig(filename=file.logs_loc+"/comparison.log", encoding='utf-8', level=logging.INFO)

class Comparison():
    def __init__(self, dnames):
        self.dnames = dnames
        self.dname = "comparison"
        self.aname = "" #a for analysis
        self.sname = ""
        self.savedir_stem = file.savedir + self.dname + "/"
        self.file_spacing = 1

    def alpha(self, fnum_range, sname = "", plot_every=100, pickle_every=None):
        self.aname = "alpha"
        self.sname = sname
        self.savedir = self.savedir_stem + self.aname
        self.pickldir = self.savedir + "/pickles"
        mkdir_if_not_exist(self.savedir)
        mkdir_if_not_exist(self.pickldir)
        if pickle_every is None:
            pickle_every = plot_every

        alpha_max = {}
        alpha_reyn = {}
        orbits = {}
        NO_VALUE = np.inf #matplot will ignore infs
        for dname in self.dnames:
            alpha_max[dname] = np.full((fnum_range[1] - fnum_range[0]), NO_VALUE)
            alpha_reyn[dname] = np.full((fnum_range[1] - fnum_range[0]), NO_VALUE)
            orbits[dname] = np.full((fnum_range[1] - fnum_range[0]), NO_VALUE)
        
            
        for i, fnum in enumerate(np.arange(fnum_range[0], fnum_range[1], self.file_spacing)):
            logging.info(f"fnum: {fnum}")
            for dname in self.dnames:
                filename = "%s/disk.out1.%05d.athdf" % (file.data_loc + dname, fnum)
                try:
                    aa = Athena_Analysis(filename=filename, grid_type=file.grid_types[dname])
                except:
                    continue

                logging.info(f"dname: {dname}")
                
                aa.get_Bfields()
                aa.get_primaries(get_rho=True, get_press=True, get_vel_r=True, get_vel_phi=True)

                stress_max = (-2/3) * aa.B_r * aa.B_phi
                stress_reyn = (-2/3) * aa.rho*aa.vel_r*aa.vel_phi
                stress_max_shell = aa.integrate(stress_max, "Shell")
                stress_reyn_shell = aa.integrate(stress_reyn, "Shell")
                pressure_shell = aa.integrate(aa.press, "Shell")
                rho_shell = aa.integrate(aa.rho, "Shell")
                alpha_max_shell = stress_max_shell / pressure_shell
                alpha_reyn_shell = stress_reyn_shell / pressure_shell

                alpha_max[dname][i] = aa.integrate(rho_shell * alpha_max_shell, "r") / aa.integrate(rho_shell, "r")
                alpha_reyn[dname][i] = aa.integrate(rho_shell * alpha_reyn_shell, "r") / aa.integrate(rho_shell, "r")

                orbits[dname][i] = aa.time / sim.binary_period

            if fnum % plot_every == 0:
                vert = 1
                horz = 2
                gs = gridspec.GridSpec(vert, horz)
                fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
                
                ax_max = fig.add_subplot(gs[0, 0])
                ax_max.set_xlabel("orbits")
                ax_max.set_ylabel(r"$\alpha$")
                ax_max.set_title(r"$\alpha_M$")
                ax_max.set_ylim([-1,1])

                ax_freyn = fig.add_subplot(gs[0, 1])
                ax_freyn.set_xlabel("orbits")
                ax_freyn.set_ylabel(r"$\alpha$")
                ax_freyn.set_title(r"$\rho v_r v_{\phi}$")

                for d, dname in enumerate(self.dnames):
                    ax_max.plot(orbits[dname], alpha_max[dname], f"C{d}-", label=dname)
                    ax_freyn.plot(orbits[dname], alpha_reyn[dname], f"C{d}-", label=dname)

                plt.legend()
                plt.subplots_adjust(top=(1-0.01*(16/vert)))
                title = r"Average $\alpha$s"
                plt.suptitle(title)
                plt.tight_layout()
                plt.savefig("%s/%s%05d%s.png" % (self.savedir, self.aname, fnum, self.sname))
                plt.close()

            if fnum % pickle_every == 0:
                with open("%s/%s%05d.dat" % (self.pickldir, self.sname, fnum), "wb") as pickle_file:
                    pickle.dump({
                        "orbits": orbits,
                        "alpha_max": alpha_max,
                        "alpha_reyn" : alpha_reyn,
                    }, pickle_file)
