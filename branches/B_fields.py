from .roots.athena_analysis import *
from .roots.misc_func import *
import pickle

logging.basicConfig(filename=file.logs_loc+"/B_fields.log", encoding='utf-8', level=logging.INFO)

class B_fields():
    def __init__(self, dname):
        self.dname = dname
        self.aname = "_B_fields" #a for analysis
        self.sname = ""
        self.data_location = file.data_loc + dname
        self.grid_type = file.grid_types[dname]
        self.is_MHD = file.MHD[dname]
        self.alpha = file.alpha[dname]
        self.savedir = file.savedir + dname + "/" + dname + self.aname
        self.pickldir = self.savedir + "/pickles"
        mkdir_if_not_exist(self.savedir)
        mkdir_if_not_exist(self.pickldir)
        self.file_spacing = find_file_spacing(self.dname)
        
    def profile(self, fnum, bar_widths=[3e0,1e1,3e-1]):
        self.sname = ""
        filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
        aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)

        aa.get_primaries(get_rho=True, get_press=True)
        aa.get_Bfields()

        rotation = -1*sim.orbital_Omega * aa.time

        vert = 2
        horz = 4
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

        ax_rho = fig.add_subplot(gs[0, 0])
        ax_rhov = fig.add_subplot(gs[1, 0])

        ax_Br = fig.add_subplot(gs[0, 1])
        ax_Brv = fig.add_subplot(gs[1, 1])

        ax_Bphi = fig.add_subplot(gs[0, 2])
        ax_Bphiv = fig.add_subplot(gs[1, 2])

        ax_Bz = fig.add_subplot(gs[0, 3])
        ax_Bzv = fig.add_subplot(gs[1, 3])

        aa.midplane_colorplot(aa.rho, ax_rho, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
        aa.midplane_colorplot(aa.rho, ax_rhov, vbound=[1e-5,1e2], slicetype='vert', rotation=rotation)
        ax_rho.set_title("Density")

        aa.midplane_colorplot(aa.B_r, ax_Br, log=False, vbound=[-bar_widths[0],bar_widths[0]], slicetype='z', rotation=rotation, cmap="PRGn")
        aa.midplane_colorplot(aa.B_r, ax_Brv, log=False, vbound=[-bar_widths[0],bar_widths[0]], slicetype='vert', rotation=rotation, cmap="PRGn")
        ax_Br.set_title(r"$B_r$")

        aa.midplane_colorplot(aa.B_phi, ax_Bphi, log=False, vbound=[-bar_widths[1],bar_widths[1]], slicetype='z', rotation=rotation, cmap="PRGn")
        aa.midplane_colorplot(aa.B_phi, ax_Bphiv, log=False, vbound=[-bar_widths[1],bar_widths[1]], slicetype='vert', rotation=rotation, cmap="PRGn")
        ax_Bphi.set_title(r"$B_{\phi}$")

        aa.midplane_colorplot(aa.B_z, ax_Bz, log=False, vbound=[-bar_widths[2],bar_widths[2]], slicetype='z', rotation=rotation, cmap="PRGn")
        aa.midplane_colorplot(aa.B_z, ax_Bzv, log=False, vbound=[-bar_widths[2],bar_widths[2]], slicetype='vert', rotation=rotation, cmap="PRGn")
        ax_Bz.set_title(r"$B_z$")

        plt.tight_layout()
        orbit = (aa.time / sim.binary_period)
        plt.subplots_adjust(top=(1-0.01*(16/vert)))
        fig.suptitle(f"orbit: {orbit:.2f}")
        plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
        plt.close()

    def vector_profile(self, fnum):
        self.sname = "vec"
        filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
        aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)

        aa.get_primaries(get_rho=True)
        aa.get_Bfields()

        rotation = -1*sim.orbital_Omega * aa.time

        B_cart = aa.native_to_cart([aa.B_z, aa.B_phi, aa.B_r])
        B_horz = B_cart[0]*np.cos(rotation) + B_cart[1]*np.sin(rotation)

        vert = 1
        horz = 2
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

        ax_rho = fig.add_subplot(gs[0, 0])
        ax_rhov = fig.add_subplot(gs[0, 1])

        aa.midplane_colorplot(aa.rho, ax_rho, vbound=[1e-5,1e2], slicetype='z')#, rotation=rotation)
        aa.midplane_colorplot(aa.rho, ax_rhov, vbound=[1e-5,1e2], slicetype='vert')#, rotation=rotation)
        ax_rho.set_title(r"Density & $\vec{B}")

        aa.midplane_vectorplot([B_cart[0], B_cart[1]], ax_rho, vbound=None, log=False, slicetype='z')#, rotation=rotation)
        aa.midplane_vectorplot([B_horz, B_cart[2]], ax_rhov, vbound=None, log=False, slicetype='vert')#, rotation=rotation)
        
        plt.tight_layout()
        orbit = (aa.time / sim.binary_period)
        plt.subplots_adjust(top=(1-0.01*(16/vert)))
        fig.suptitle(f"{self.dname} orbit: {orbit:.2f}")
        plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
        plt.close()

    def pressure(self, fnum):
        self.sname = "pressure"
        filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
        aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)

        aa.get_primaries(get_rho=True, get_press=True)
        aa.get_Bfields()

        B_press = (aa.B_r*aa.B_r + aa.B_phi*aa.B_phi + aa.B_z*aa.B_z)/2

        rotation = -1*sim.orbital_Omega * aa.time

        vert = 2
        horz = 3
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

        ax_rho = fig.add_subplot(gs[0, 0])
        ax_press = fig.add_subplot(gs[1, 0])

        ax_Bpress = fig.add_subplot(gs[0, 1])
        ax_Br = fig.add_subplot(gs[1, 1])

        ax_Bphi = fig.add_subplot(gs[0, 2])
        ax_Bz = fig.add_subplot(gs[1, 2])

        aa.midplane_colorplot(aa.rho, ax_rho, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
        aa.midplane_colorplot(aa.press, ax_press, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
        ax_rho.set_title("Density")
        ax_press.set_title("Presure")

        aa.midplane_colorplot(B_press, ax_Bpress, log=True, vbound=[1e-5, 1e2], slicetype='z', rotation=rotation)
        aa.midplane_colorplot((aa.B_r*aa.B_r)/2, ax_Br, log=True, vbound=[1e-5, 1e2], slicetype='z', rotation=rotation)
        aa.midplane_colorplot((aa.B_phi*aa.B_phi)/2, ax_Bphi, log=True, vbound=[1e-5, 1e2], slicetype='z', rotation=rotation)
        aa.midplane_colorplot((aa.B_z*aa.B_z)/2, ax_Bz, log=True, vbound=[1e-5, 1e2], slicetype='z', rotation=rotation)
        ax_Bpress.set_title("Magnetic Pressure")
        ax_Br.set_title(r"$\frac{1}{2} B_{r}^2$")
        ax_Bphi.set_title(r"$\frac{1}{2} B_{\phi}^2$")
        ax_Bz.set_title(r"$\frac{1}{2} B_{z}^2$")

        plt.tight_layout()
        orbit = (aa.time / sim.binary_period)
        plt.subplots_adjust(top=(1-0.01*(16/vert)))
        fig.suptitle(f"orbit: {orbit:.2f}")
        plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
        plt.close()

    def alfven(self, fnum):
        self.sname = "alfven"
        filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
        aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)

        aa.native_grid(get_r=True, get_z=True)
        aa.get_primaries(get_rho=True, get_press=True)
        aa.get_Bfields()

        B_press = (aa.B_r*aa.B_r + aa.B_phi*aa.B_phi + aa.B_z*aa.B_z)/2

        alfven_speed = np.sqrt(B_press/aa.rho)
        sound_speed = np.sqrt((5/3)*(aa.press/aa.rho))
        height = aa.possible_z[-1] - aa.possible_z[0]
        local_omega = np.sqrt(sim.gm1/(aa.r*aa.r*aa.r))

        rotation = -1*sim.orbital_Omega * aa.time

        vert = 2
        horz = 2
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

        ax_rho = fig.add_subplot(gs[0, 0])
        ax_alfven = fig.add_subplot(gs[1, 0])
        ax_alfven_sound = fig.add_subplot(gs[1, 1])
        ax_alfven_height = fig.add_subplot(gs[0, 1])

        aa.midplane_colorplot(aa.rho, ax_rho, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
        aa.midplane_colorplot(alfven_speed, ax_alfven, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
        divnorm = colors.TwoSlopeNorm(vmin=0, vcenter=1, vmax=2)
        aa.midplane_colorplot(alfven_speed/sound_speed, ax_alfven_sound, slicetype='z', norm=divnorm, cmap="PRGn", rotation=rotation)
        aa.midplane_colorplot((2*np.pi/np.sqrt(3))*alfven_speed/(local_omega*height), ax_alfven_height, slicetype='z', cmap="PRGn", norm=divnorm, rotation=rotation)
        ax_rho.set_title("Density")
        ax_alfven.set_title("Alfven Speed")
        ax_alfven_sound.set_title(r"$\frac{\Omega H}{v_{S}}$")
        ax_alfven_height.set_title(r"$\frac{v_{A}}{\Omega H} \frac{2\pi}{\sqrt{3}}$")

        plt.tight_layout()
        orbit = (aa.time / sim.binary_period)
        plt.subplots_adjust(top=(1-0.01*(16/vert)))
        fig.suptitle(f"orbit: {orbit:.2f}")
        plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
        plt.close()
    
    def beta(self, fnum):
        self.sname = "beta"
        filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
        aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)

        aa.get_primaries(get_rho=True, get_press=True)
        aa.get_Bfields()

        B_press = (aa.B_r*aa.B_r + aa.B_phi*aa.B_phi + aa.B_z*aa.B_z)/2

        rotation = -1*sim.orbital_Omega * aa.time

        vert = 2
        horz = 2
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

        ax_rho = fig.add_subplot(gs[0, 0])
        ax_press = fig.add_subplot(gs[1, 0])

        ax_Bpress = fig.add_subplot(gs[0, 1])
        ax_beta = fig.add_subplot(gs[1, 1])

        aa.midplane_colorplot(aa.rho, ax_rho, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
        aa.midplane_colorplot(aa.press, ax_press, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
        ax_rho.set_title("Density")
        ax_press.set_title("Pressure")

        aa.midplane_colorplot(B_press, ax_Bpress, log=True, vbound=[1e-5, 1e2], slicetype='z', rotation=rotation)
        aa.midplane_colorplot(aa.press / B_press, ax_beta, log=True, vbound=[1e-5, 1e3], slicetype='z', rotation=rotation)
        ax_Bpress.set_title("Magnetic Pressure")
        ax_beta.set_title(r"$\beta$")

        plt.tight_layout()
        orbit = (aa.time / sim.binary_period)
        plt.subplots_adjust(top=(1-0.01*(16/vert)))
        fig.suptitle(f"{self.dname} orbit: {orbit:.2f}")
        plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
        plt.close()

    def alpha_plot(self, fnum):
        self.sname = "alpha"
        filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
        aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)

        aa.get_primaries(get_rho=True, get_press=True)
        aa.get_Bfields()

        B_press = (aa.B_r*aa.B_r + aa.B_phi*aa.B_phi + aa.B_z*aa.B_z)/2

        alpha = (- aa.B_r * aa.B_phi) / aa.press

        rotation = -1*sim.orbital_Omega * aa.time

        vert = 2
        horz = 2
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

        ax_rho = fig.add_subplot(gs[0, 0])
        ax_press = fig.add_subplot(gs[1, 0])

        ax_alpha = fig.add_subplot(gs[0, 1])
        ax_beta = fig.add_subplot(gs[1, 1])

        aa.midplane_colorplot(aa.rho, ax_rho, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
        aa.midplane_colorplot(B_press, ax_press, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
        ax_rho.set_title("Density")
        ax_press.set_title("Magnetic Pressure")

        aa.midplane_colorplot(alpha, ax_alpha, log=True, vbound=[-1e2, 1e2], slicetype='z', rotation=rotation)
        aa.midplane_colorplot(aa.press / B_press, ax_beta, log=True, vbound=[1e-5, 1e3], slicetype='z', rotation=rotation)
        ax_alpha.set_title(r"$\alpha_M$")
        ax_beta.set_title(r"$\beta$")

        plt.tight_layout()
        orbit = (aa.time / sim.binary_period)
        plt.subplots_adjust(top=(1-0.01*(16/vert)))
        fig.suptitle(f"{self.dname} orbit: {orbit:.2f}")
        plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
        plt.close()

    def evolution(self, fnum_range, plot_every=100, pickle_every=None):
        self.sname = "evolution"
        fnum_range = np.arange(fnum_range[0], fnum_range[1], self.file_spacing)
        if pickle_every is None:
            pickle_every = plot_every

        B_fields = np.zeros([3, len(fnum_range)])
        B_pressures = np.zeros([3, len(fnum_range)])
        orbits = np.zeros(len(fnum_range))

        for fnum, i in enumerate(fnum_range):
            filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
            try:
                aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)
            except:
                vert = 2
                horz = 3
                gs = gridspec.GridSpec(vert, horz)
                fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

                ax_Br = fig.add_subplot(gs[0, 0])
                ax_Brp = fig.add_subplot(gs[1, 0])
                ax_Br.set_title(r"$B_r$")
                ax_Brp.set_title(r"$\frac{1}{2} B_r^2$")

                ax_Bp = fig.add_subplot(gs[0, 1])
                ax_Bpp = fig.add_subplot(gs[1, 1])
                ax_Bp.set_title(r"$B_{\phi}$")
                ax_Bpp.set_title(r"$\frac{1}{2} B_{\phi}^2$")

                ax_Bz = fig.add_subplot(gs[0, 2])
                ax_Bzp = fig.add_subplot(gs[1, 2])
                ax_Bz.set_title(r"$B_z$")
                ax_Bzp.set_title(r"$\frac{1}{2} B_z^2$")

                ax_Br.plot(orbits[:i], B_fields[0,:i], c="b")
                ax_Bp.plot(orbits[:i], B_fields[1,:i], c="b")
                ax_Bz.plot(orbits[:i], B_fields[2,:i], c="b")

                ax_Brp.plot(orbits[:i], B_pressures[0,:i], c="m")
                ax_Bpp.plot(orbits[:i], B_pressures[1,:i], c="m")
                ax_Bzp.plot(orbits[:i], B_pressures[2,:i], c="m")

                plt.tight_layout()
                
                plt.subplots_adjust(top=(1-0.01*(16/vert)))
                fig.suptitle(f"{self.dname}")
                plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
                plt.close()
                return

            aa.get_Bfields()
            
            B_fields[0, i] = aa.integrate(aa.B_r, "All")
            B_fields[1, i] = aa.integrate(aa.B_phi, "All")
            B_fields[2, i] = aa.integrate(aa.B_z, "All")

            B_pressures[0, i] = aa.integrate(0.5*aa.B_r*aa.B_r, "All")
            B_pressures[1, i] = aa.integrate(0.5*aa.B_phi*aa.B_phi, "All")
            B_pressures[2, i] = aa.integrate(0.5*aa.B_z*aa.B_z, "All")

            orbits[i] = (aa.time / sim.binary_period)

            if (i % plot_every == 0):
                vert = 2
                horz = 3
                gs = gridspec.GridSpec(vert, horz)
                fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

                ax_Br = fig.add_subplot(gs[0, 0])
                ax_Brp = fig.add_subplot(gs[1, 0])
                ax_Br.set_title(r"$B_r$")
                ax_Brp.set_title(r"$\frac{1}{2} B_r^2$")

                ax_Bp = fig.add_subplot(gs[0, 1])
                ax_Bpp = fig.add_subplot(gs[1, 1])
                ax_Bp.set_title(r"$B_{\phi}$")
                ax_Bpp.set_title(r"$\frac{1}{2} B_{\phi}^2$")

                ax_Bz = fig.add_subplot(gs[0, 2])
                ax_Bzp = fig.add_subplot(gs[1, 2])
                ax_Bz.set_title(r"$B_z$")
                ax_Bzp.set_title(r"$\frac{1}{2} B_z^2$")

                ax_Br.plot(orbits[:i+1], B_fields[0,:i+1], c="b")
                ax_Bp.plot(orbits[:i+1], B_fields[1,:i+1], c="b")
                ax_Bz.plot(orbits[:i+1], B_fields[2,:i+1], c="b")

                ax_Brp.plot(orbits[:i+1], B_pressures[0,:i+1], c="m")
                ax_Bpp.plot(orbits[:i+1], B_pressures[1,:i+1], c="m")
                ax_Bzp.plot(orbits[:i+1], B_pressures[2,:i+1], c="m")

                plt.tight_layout()
                
                plt.subplots_adjust(top=(1-0.01*(16/vert)))
                fig.suptitle(f"{self.dname}")
                plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
                plt.close()

            if i % pickle_every == 0:
                with open("%s/%s%05d.dat" % (self.pickldir, self.sname, fnum), "wb") as pickle_file:
                    pickle.dump({
                        "B_fields": B_fields[:,:i+1],
                        "B_pressures": B_pressures[:,:i+1],
                        "orbits": orbits[:i+1]
                    }, pickle_file)

    def alpha_beta_evolution(self, fnum_range, plot_every=100, pickle_every=None):
        self.sname = "alpha_beta_evol"
        fnum_range = np.arange(fnum_range[0], fnum_range[1], self.file_spacing)
        if pickle_every is None:
            pickle_every = plot_every

        alpha = np.zeros([len(fnum_range)])
        beta = np.zeros([len(fnum_range)])
        rho = np.zeros([len(fnum_range)])
        press = np.zeros([len(fnum_range)])
        orbits = np.zeros(len(fnum_range))

        for fnum, i in enumerate(fnum_range):
            filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
            try:
                aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)
            except:
                vert = 2
                horz = 2
                gs = gridspec.GridSpec(vert, horz)
                fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

                ax_alpha = fig.add_subplot(gs[0, 0])
                ax_beta = fig.add_subplot(gs[0, 1])
                ax_rho = fig.add_subplot(gs[1, 0])
                ax_press = fig.add_subplot(gs[1, 1])
                ax_alpha.set_title(r"Ave $\alpha$")
                ax_beta.set_title(r"Ave $\beta$")
                ax_rho.set_title(r"$\rho$")
                ax_press.set_title("Ave P")

                ax_alpha.plot(orbits[:i], alpha[:i], c="m")
                ax_beta.plot(orbits[:i], beta[:i], c="m")
                ax_rho.plot(orbits[:i], rho[:i], c="b")
                ax_press.plot(orbits[:i], press[:i], c="b")
                
                plt.tight_layout()
                
                plt.subplots_adjust(top=(1-0.01*(16/vert)))
                fig.suptitle(f"{self.dname}")
                plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
                plt.close()
                return

            aa.get_Bfields()
            aa.get_primaries(get_rho=True, get_press=True)
            
            volume = aa.integrate(1, "All")
            rho[i] = aa.integrate(aa.rho, "All")
            press[i] = aa.integrate(aa.press, "All") / volume
            alpha[i] = (aa.integrate(aa.B_r * aa.B_phi, "All") / volume) / press[i]
            beta[i] = press[i] / (aa.integrate(0.5*(aa.B_r*aa.B_r + aa.B_phi*aa.B_phi + aa.B_z*aa.B_z), "All") / volume)
            orbits[i] = (aa.time / sim.binary_period)

            if (i % plot_every == 0):
                vert = 2
                horz = 2
                gs = gridspec.GridSpec(vert, horz)
                fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

                ax_alpha = fig.add_subplot(gs[0, 0])
                ax_beta = fig.add_subplot(gs[0, 1])
                ax_rho = fig.add_subplot(gs[1, 0])
                ax_press = fig.add_subplot(gs[1, 1])
                ax_alpha.set_title(r"Ave $\alpha$")
                ax_beta.set_title(r"Ave $\beta$")
                ax_rho.set_title(r"$\rho$")
                ax_press.set_title("Ave P")

                ax_alpha.plot(orbits[:i], alpha[:i], c="m")
                ax_beta.plot(orbits[:i], beta[:i], c="m")
                ax_rho.plot(orbits[:i], rho[:i], c="b")
                ax_press.plot(orbits[:i], press[:i], c="b")
                
                plt.tight_layout()
                
                plt.subplots_adjust(top=(1-0.01*(16/vert)))
                fig.suptitle(f"{self.dname}")
                plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
                plt.close()

            if i % pickle_every == 0:
                with open("%s/%s%05d.dat" % (self.pickldir, self.sname, fnum), "wb") as pickle_file:
                    pickle.dump({
                        "alpha": alpha[:,:i+1],
                        "beta": beta[:,:i+1],
                        "orbits": orbits[:i+1]
                    }, pickle_file)

    def alpha_beta_radial(self, fnum):
        self.sname = "alpha_beta_rad"

        filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
        aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)

        aa.get_Bfields()
        aa.get_primaries(get_rho=True, get_press=True)
        aa.native_grid()
        
        shells = aa.integrate(1, "Shell")
        rho = aa.integrate(aa.rho, "Shell") / shells
        press = aa.integrate(aa.press, "Shell") / shells
        alpha = (aa.integrate(-1 *aa.B_r * aa.B_phi, "Shell") / shells) / press
        beta = press / (aa.integrate(0.5*(aa.B_r*aa.B_r + aa.B_phi*aa.B_phi + aa.B_z*aa.B_z), "Shell") / shells)
        r_axis = aa.possible_r

        vert = 2
        horz = 2
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

        ax_alpha = fig.add_subplot(gs[0, 0])
        ax_beta = fig.add_subplot(gs[0, 1])
        ax_rho = fig.add_subplot(gs[1, 0])
        ax_press = fig.add_subplot(gs[1, 1])
        ax_alpha.set_title(r"Ave $\alpha$")
        ax_beta.set_title(r"Ave $\beta$")
        ax_rho.set_title(r"Ave $\rho$")
        ax_press.set_title("Ave P")
        ax_beta.set_yscale("log")
        ax_rho.set_yscale("log")
        ax_press.set_yscale("log")

        ax_alpha.plot(r_axis, alpha, c="m")
        ax_beta.plot(r_axis, beta, c="m")
        ax_rho.plot(r_axis, rho, c="b")
        ax_press.plot(r_axis, press, c="b")
        
        plt.tight_layout()
        
        plt.subplots_adjust(top=(1-0.01*(16/vert)))
        fig.suptitle(f"{self.dname}")
        plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
        plt.close()