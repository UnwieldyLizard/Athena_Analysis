from .roots.athena_analysis import *
from .roots.misc_func import *
import pickle

logging.basicConfig(filename=file.logs_loc+"/comparison.log", encoding='utf-8', level=logging.INFO)

class Comparison():
    def __init__(self, dnames, file_spacing):
        self.dnames = dnames
        self.dname = "comparison"
        self.aname = "" #a for analysis
        self.sname = ""
        self.savedir_stem = file.savedir + self.dname + "/"
        self.file_spacing = file_spacing

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

                alpha_max[dname][i] = np.sum(rho_shell * alpha_max_shell) / np.sum(rho_shell)
                alpha_reyn[dname][i] = np.sum(rho_shell * alpha_reyn_shell) / np.sum(rho_shell)

                orbits[dname][i] = aa.time / sim.binary_period

            if fnum % plot_every == 0:
                vert = 1
                horz = 2
                gs = gridspec.GridSpec(vert, horz)
                fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
                
                ax_max = fig.add_subplot(gs[0, 0])
                ax_max.set_xlabel("orbits")
                ax_max.set_ylabel(r"$\alpha$")
                ax_max.set_title(r"$\langle\alpha_M\rangle$")
                ax_max.set_ylim([-1,1])

                ax_freyn = fig.add_subplot(gs[0, 1])
                ax_freyn.set_xlabel("orbits")
                ax_freyn.set_ylabel(r"$\alpha$")
                ax_freyn.set_title(r"$\langle\alpha_R\rangle$")

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
    
    def alpha_replot(self, file, sname = "", log=True, ylims=[None, None], cutoffmin=35):
        self.aname = "alpha"
        self.sname = sname
        self.savedir = self.savedir_stem + self.aname
        self.pickldir = self.savedir + "/pickles"

        orbit_ave_width = 1

        with open("%s/%s%s.dat" % (self.pickldir, self.sname, file), "rb") as pickle_file:
            data = pickle.load(pickle_file)

        vert = 1
        horz = 2
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
        
        ax_m = fig.add_subplot(gs[0, 0])
        ax_m.set_xlabel("orbits")
        ax_m.set_ylabel(r"$\alpha$")
        ax_m.set_title(r"$\langle\alpha_M\rangle$")
        if log:
            ax_m.set_yscale("log")
        if ylims[0] is not None:
            ax_m.set_ylim(ylims[0])

        for d, dname in enumerate(self.dnames):
            box_ave = self.box_average(data["alpha_max"][dname], orbit_ave_width * sim.filenums_per_orbit)

            arg = np.argmin(abs(data["orbits"][dname] - cutoffmin))
            filtered_alpha = np.ma.compressed(np.ma.masked_where(data["alpha_max"][dname][arg:] == np.inf, data["alpha_max"][dname][arg:]))
            ave = np.sum(filtered_alpha) / len(filtered_alpha)
            ax_m.plot(data["orbits"][dname], box_ave, f"C{d}-", label=dname+f": {ave:.3f}", linewidth=1)
            ax_m.plot(data["orbits"][dname], np.full(data["orbits"][dname].shape ,ave), f"C{d}--", linewidth=1)

        ax_m.legend()
        
        '''
        ax_r = fig.add_subplot(gs[0, 1])
        ax_r.set_xlabel("orbits")
        ax_r.set_ylabel(r"$\alpha$")
        ax_r.set_title(r"$\rho v_r v_{\phi}$")
        if log:
            ax_r.set_yscale("log")
        if ylims[1] is not None:
            ax_r.set_ylim(ylims[1])

        for d, dname in enumerate(self.dnames):
            arg = np.argmin(abs(data["orbits"][dname] - cutoffmin))
            filtered_alpha = np.ma.compressed(np.ma.masked_where(data["alpha_reyn"][dname][arg:] == np.inf, data["alpha_reyn"][dname][arg:]))
            ave = np.sum(filtered_alpha) / len(filtered_alpha)
            ax_r.plot(data["orbits"][dname], data["alpha_reyn"][dname], f"C{d}-", label=dname+f": {ave:.3f}", linewidth=0.5)
            ax_r.plot(data["orbits"][dname], np.full(data["orbits"][dname].shape ,ave), f"C{d}--", linewidth=1)
        
        ax_r.legend()
        '''

        ax_rb = fig.add_subplot(gs[0, 1])
        ax_rb.set_xlabel("orbits")
        ax_rb.set_ylabel(r"$\alpha$")
        ax_rb.set_title(r"$\langle\alpha_R\rangle$")
        if log:
            ax_rb.set_yscale("log")
        if ylims[2] is not None:
            ax_rb.set_ylim(ylims[2])

        for d, dname in enumerate(self.dnames):
            box_ave = self.box_average(data["alpha_reyn"][dname], orbit_ave_width * sim.filenums_per_orbit)
            
            arg = np.argmin(abs(data["orbits"][dname] - cutoffmin))
            filtered_alpha = np.ma.compressed(np.ma.masked_where(data["alpha_reyn"][dname][arg:] == np.inf, data["alpha_reyn"][dname][arg:]))
            ave = np.sum(filtered_alpha) / len(filtered_alpha)
            ax_rb.plot(data["orbits"][dname], box_ave, f"C{d}-", label=dname+f": {ave:.3f}", linewidth=0.5)
            ax_rb.plot(data["orbits"][dname], np.full(data["orbits"][dname].shape ,ave), f"C{d}--", linewidth=1)
        
        ax_rb.legend()

        plt.tight_layout()
        plt.savefig("%s/%s%s%s_replot.png" % (self.savedir, self.aname, file, self.sname))
        plt.close()

    def box_average(self, data, filenum_width):
        box_ave = np.copy(data)
        for i in range(len(data)):
            if i-round(0.5*filenum_width) < 0:
                box_ave[i] = np.sum(data[:i+round(0.5*filenum_width)]) / (i+round(0.5*filenum_width))
                continue
            if i+round(0.5*filenum_width) > len(data) - 1:
                box_ave[i] = np.sum(data[i-round(0.5*filenum_width):]) / (len(data) - (i-round(0.5*filenum_width)))
                continue
            box_ave[i] = np.sum(data[i-round(0.5*filenum_width):i+round(0.5*filenum_width)]) / filenum_width
        
        return box_ave

    def beta(self, fnum_range, sname = "", plot_every=100, pickle_every=None):
        self.aname = "beta"
        self.sname = sname
        self.savedir = self.savedir_stem + self.aname
        self.pickldir = self.savedir + "/pickles"
        mkdir_if_not_exist(self.savedir)
        mkdir_if_not_exist(self.pickldir)
        if pickle_every is None:
            pickle_every = plot_every

        beta = {}
        orbits = {}

        NO_VALUE = np.inf #matplot will ignore infs
        for dname in self.dnames:
            beta[dname] = np.full((fnum_range[1] - fnum_range[0]), NO_VALUE)
            orbits[dname] = np.full((fnum_range[1] - fnum_range[0]), NO_VALUE)

        finished_dnames = []
            
        for i, fnum in enumerate(np.arange(fnum_range[0], fnum_range[1], self.file_spacing)):
            logging.info(f"fnum: {fnum}")
            for dname in self.dnames:
                if dname in finished_dnames:
                    continue

                filename = "%s/disk.out1.%05d.athdf" % (file.data_loc + dname, fnum)
                try:
                    aa = Athena_Analysis(filename=filename, grid_type=file.grid_types[dname])
                except:
                    finished_dnames.append(dname)
                    if len(finished_dnames) == len(self.dnames):
                        break
                    continue

                logging.info(f"dname: {dname}")
                
                aa.get_Bfields()
                aa.get_primaries(get_press=True, get_rho=True)

                B_press = (aa.B_r*aa.B_r + aa.B_phi*aa.B_phi + aa.B_z*aa.B_z)/2
                ave_B_press = aa.integrate(B_press, "All")
                ave_press = aa.integrate(aa.press, "All")
                
                if ave_B_press == 0:
                    beta[dname][i] = np.inf
                else:
                    beta[dname][i] = ave_press / ave_B_press

                orbits[dname][i] = aa.time / sim.binary_period

            if fnum % plot_every == 0:
                vert = 1
                horz = 1
                gs = gridspec.GridSpec(vert, horz)
                fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
                
                ax = fig.add_subplot(gs[0, 0])
                ax.set_xlabel("orbits")
                ax.set_ylabel(r"$\beta$")
                ax.set_title(r"$\langle\beta\rangle$")

                for d, dname in enumerate(self.dnames):
                    ax.plot(orbits[dname], beta[dname], f"C{d}-", label=dname)
                
                plt.legend()
                plt.tight_layout()
                plt.savefig("%s/%s%05d%s.png" % (self.savedir, self.aname, fnum, self.sname))
                plt.close()

            if fnum % pickle_every == 0:
                with open("%s/%s%05d.dat" % (self.pickldir, self.sname, fnum), "wb") as pickle_file:
                    pickle.dump({
                        "orbits": orbits,
                        "beta": beta,
                    }, pickle_file)

        #final plot
        vert = 1
        horz = 1
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
        
        ax = fig.add_subplot(gs[0, 0])
        ax.set_xlabel("orbits")
        ax.set_ylabel(r"$\beta$")
        ax.set_title(r"Average $\beta$")

        for d, dname in enumerate(self.dnames):
            ax.plot(orbits[dname], beta[dname], f"C{d}-", label=dname)
        
        plt.legend()
        plt.tight_layout()
        plt.savefig("%s/%s%s%s.png" % (self.savedir, self.aname, "full", self.sname))
        plt.close()

        #final pickl
        with open("%s/%s%s.dat" % (self.pickldir, self.sname, "full"), "wb") as pickle_file:
            pickle.dump({
                "orbits": orbits,
                "beta": beta,
            }, pickle_file)

    def beta_replot(self, file, sname = "", log=True, ylims=None, cutoffmin=35):
        self.aname = "beta"
        self.sname = sname
        self.savedir = self.savedir_stem + self.aname
        self.pickldir = self.savedir + "/pickles"

        orbit_ave_width = 1

        with open("%s/%s%s.dat" % (self.pickldir, self.sname, file), "rb") as pickle_file:
            data = pickle.load(pickle_file)

        vert = 1
        horz = 1
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
        
        ax = fig.add_subplot(gs[0, 0])
        ax.set_xlabel("orbits")
        ax.set_ylabel(r"$\beta$")
        ax.set_title(r"$\langle\beta\rangle$")
        if log:
            ax.set_yscale("log")
        if ylims is not None:
            ax.set_ylim(ylims)

        for d, dname in enumerate(self.dnames):
            box_ave = self.box_average(data["beta"][dname], orbit_ave_width * sim.filenums_per_orbit)

            arg = np.argmin(abs(data["orbits"][dname] - cutoffmin))
            filtered_beta = np.ma.compressed(np.ma.masked_where(data["beta"][dname][arg:] == np.inf, data["beta"][dname][arg:]))
            ave = np.sum(filtered_beta) / len(filtered_beta)
            ax.plot(data["orbits"][dname], box_ave, f"C{d}-", label=dname+f": {ave:.2f}", linewidth=1)
            ax.plot(data["orbits"][dname], np.full(data["orbits"][dname].shape ,ave), f"C{d}--", linewidth=1)
        
        plt.legend()
        plt.tight_layout()
        plt.savefig("%s/%s%s%s_replot.png" % (self.savedir, self.aname, file, self.sname))
        plt.close()