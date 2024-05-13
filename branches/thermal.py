from .roots.athena_analysis import *
from .roots.misc_func import *
from datetime import datetime
import gc
import pickle
#from pympler.tracker import SummaryTracker
#tracker = SummaryTracker()

logging.basicConfig(filename=file.logs_loc+"/momentum.log", encoding='utf-8', level=logging.INFO)

class Thermal():
    def __init__(self, dname):
        self.dname = dname
        self.aname = "_thermal" #a for analysis
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

    def entropy_profile(self, fnum):
        now = datetime.now()
        self.sname = "entropy_profile"

        logging.info("fnum = %d" % fnum)
        filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
        aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)

        aa.get_primaries(get_rho=True, get_press=True)

        #Transforming to inertial frame
        rotation = -1*sim.orbital_Omega * aa.time

        vert = 2
        horz = 4
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

        #'''
        ax_rho = fig.add_subplot(gs[0, 0])
        ax_rhov = fig.add_subplot(gs[1, 0])

        ax_press = fig.add_subplot(gs[0, 1])
        ax_pressv = fig.add_subplot(gs[1, 1])

        ax_3 = fig.add_subplot(gs[0, 2])
        ax_3v = fig.add_subplot(gs[1, 2])

        ax_4 = fig.add_subplot(gs[0, 3])
        ax_4v = fig.add_subplot(gs[1, 3])

        aa.midplane_colorplot(aa.rho, ax_rho, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
        aa.midplane_colorplot(aa.rho, ax_rhov, vbound=[1e-5,1e2], slicetype='vert', rotation=rotation)
        ax_rho.set_title("Density")

        aa.midplane_colorplot(aa.press, ax_press, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
        aa.midplane_colorplot(aa.press, ax_pressv, vbound=[1e-5,1e2], slicetype='vert', rotation=rotation)
        ax_press.set_title("Pressure")

        aa.midplane_colorplot(aa.press / aa.rho, ax_3, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
        aa.midplane_colorplot(aa.press / aa.rho, ax_3v, vbound=[1e-5,1e2], slicetype='vert', rotation=rotation)
        ax_3.set_title(r"$\frac{P}{\rho}$")

        gamma = 5/3

        aa.midplane_colorplot(aa.press / (aa.rho ** gamma), ax_4, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
        aa.midplane_colorplot(aa.press / (aa.rho ** gamma), ax_4v, vbound=[1e-5,1e2], slicetype='vert', rotation=rotation)
        ax_4.set_title(r"$\frac{P}{\rho^{\frac{5}{3}}}$")

        plt.tight_layout()
        orbit = (aa.time / sim.binary_period)
        plt.subplots_adjust(top=(1-0.01*(16/vert)))
        fig.suptitle(f"orbit: {orbit:.2f}")
        plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
        plt.close()
        del aa
        gc.collect()
        logging.info(datetime.now()-now)
    
    def evolution(self, fnum_range, plot_every=100, pickle_every=None):
        self.sname = "evolution"
        fnum_range = np.arange(fnum_range[0], fnum_range[1], self.file_spacing)
        if pickle_every is None:
            pickle_every = plot_every

        momenta = np.zeros([3, len(fnum_range)])
        energies = np.zeros([3, len(fnum_range)])
        orbits = np.zeros(len(fnum_range))

        for fnum, i in enumerate(fnum_range):
            filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
            aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)

            aa.native_grid(get_r=True)
            aa.get_primaries(get_rho=True, get_vel_r=True, get_vel_phi=True, get_vel_z=True)
            
            v_phi = aa.vel_phi - (np.sqrt(sim.gm1/aa.r)) + aa.r*sim.orbital_Omega
            
            momenta[0, i] = aa.integrate(aa.rho * aa.vel_r, "All")
            momenta[1, i] = aa.integrate(aa.rho * v_phi, "All")
            momenta[2, i] = aa.integrate(aa.rho * aa.vel_z, "All")

            energies[0, i] = aa.integrate(0.5*aa.vel_r*aa.vel_r, "All")
            energies[1, i] = aa.integrate(0.5*v_phi*v_phi, "All")
            energies[2, i] = aa.integrate(0.5*aa.vel_z*aa.vel_z, "All")

            orbits[i] = (aa.time / sim.binary_period)

            if i % plot_every == 0:
                vert = 2
                horz = 3
                gs = gridspec.GridSpec(vert, horz)
                fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

                ax_Br = fig.add_subplot(gs[0, 0])
                ax_Brp = fig.add_subplot(gs[1, 0])
                ax_Br.set_title(r"$m v_r$")
                ax_Brp.set_title(r"$\frac{1}{2} m v_r^2$")

                ax_Bp = fig.add_subplot(gs[0, 1])
                ax_Bpp = fig.add_subplot(gs[1, 1])
                ax_Bp.set_title(r"$m v_{\phi}$")
                ax_Bpp.set_title(r"$\frac{1}{2} m v_{\phi}^2$")

                ax_Bz = fig.add_subplot(gs[0, 2])
                ax_Bzp = fig.add_subplot(gs[1, 2])
                ax_Bz.set_title(r"$m v_z$")
                ax_Bzp.set_title(r"$\frac{1}{2} m v_z^2$")

                ax_Br.plot(orbits[:i+1], momenta[0,:i+1], c="b")
                ax_Bp.plot(orbits[:i+1], momenta[1,:i+1], c="b")
                ax_Bz.plot(orbits[:i+1], momenta[2,:i+1], c="b")

                ax_Brp.plot(orbits[:i+1], energies[0,:i+1], c="m")
                ax_Bpp.plot(orbits[:i+1], energies[1,:i+1], c="m")
                ax_Bzp.plot(orbits[:i+1], energies[2,:i+1], c="m")

                plt.tight_layout()
                
                plt.subplots_adjust(top=(1-0.01*(16/vert)))
                fig.suptitle(f"{self.dname}")
                plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
                plt.close()

            if i % pickle_every == 0:
                with open("%s/%s%05d.dat" % (self.pickldir, self.sname, fnum), "wb") as pickle_file:
                    pickle.dump({
                        "momenta": momenta[:,:i+1],
                        "energies": energies[:,:i+1],
                        "orbits": orbits[:i+1]
                    }, pickle_file)
    
    def radial_profile(self, fnum):
        self.sname = "radial"

        filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
        aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)

        aa.get_primaries(get_rho=True, get_vel_r=True, get_vel_phi=True, get_vel_z=True)
        aa.native_grid(get_r=True)
        
        v_phi = aa.vel_phi - (np.sqrt(sim.gm1/aa.r)) + aa.r*sim.orbital_Omega
        
        shells = aa.integrate(1, "Shell")
        
        r_momentum = aa.integrate(aa.vel_r * aa.rho, "Shell") / shells
        phi_momentum = aa.integrate(v_phi * aa.rho, "Shell") / shells
        z_momentum = aa.integrate(aa.vel_z * aa.rho, "Shell") / shells
        r_energy = aa.integrate(0.5 * aa.vel_r * aa.vel_r * aa.rho, "Shell") / shells
        phi_energy = aa.integrate(0.5 * v_phi * v_phi * aa.rho, "Shell") / shells
        z_energy = aa.integrate(0.5 * aa.vel_z * aa.vel_z * aa.rho, "Shell") / shells
        r_axis = aa.possible_r

        vert = 2
        horz = 3
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

        ax_mr = fig.add_subplot(gs[0, 0])
        ax_mphi = fig.add_subplot(gs[0, 1])
        ax_mz = fig.add_subplot(gs[0, 2])
        
        ax_er = fig.add_subplot(gs[1, 0])
        ax_ephi = fig.add_subplot(gs[1, 1])
        ax_ez = fig.add_subplot(gs[1, 2])
        
        ax_mr.set_title(r"Ave $m v_r$")
        ax_mphi.set_title(r"Ave $m v_{\phi}$")
        ax_mz.set_title(r"Ave $m v_z$")
        ax_er.set_title(r"Ave $\frac{1}{2} m v_r^2$")
        ax_ephi.set_title(r"Ave $\frac{1}{2} m v_{\phi}^2$")
        ax_ez.set_title(r"Ave $\frac{1}{2} m v_z^2$")
        
        ax_mr.plot(r_axis, r_momentum, c="b")
        ax_mphi.plot(r_axis, phi_momentum, c="b")
        ax_mz.plot(r_axis, z_momentum, c="b")
        ax_er.plot(r_axis, r_energy, c="m")
        ax_ephi.plot(r_axis, phi_energy, c="m")
        ax_ez.plot(r_axis, z_energy, c="m")

        plt.tight_layout()

        orbit = (aa.time / sim.binary_period)
        
        plt.subplots_adjust(top=(1-0.01*(16/vert)))
        fig.suptitle(f"{self.dname} orbit: {orbit:.2f}")
        plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
        plt.close()

    def temperature_profile(self, fnum):
        self.sname = "temp"

        filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)
        aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)

        aa.get_primaries(get_rho=True, get_press=True)
        temp = (aa.press/aa.rho) * (1) #put R/m here is I figure it out

        rad_temp = aa.integrate(temp, "shell") / aa.integrate(1, "shell")

        vert = 1
        horz = 2
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

        ax_temp = fig.add_subplot(gs[0, 0])
        ax_rtemp = fig.add_subplot(gs[0, 1])
        
        aa.midplane_colorplot(temp, ax_temp)
        ax_temp.set_title(r"$\left(\frac{R}{m}\right)T$")

        ax_rtemp.plot(aa.possible_r, rad_temp)
        ax_rtemp.plot(aa.possible_r, 32.5/aa.possible_r, linestyle="--", label=r"$\frac{32.5}{r}$")
        ax_rtemp.set_title(r"$\left(\frac{R}{m}\right)T$")
        ax_rtemp.set_xlabel("r")
        ax_rtemp.legend()

        plt.tight_layout()

        orbit = (aa.time / sim.binary_period)
        
        plt.subplots_adjust(top=(1-0.01*(16/vert)))
        fig.suptitle(f"{self.dname} orbit: {orbit:.2f}")
        plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, self.sname))
        plt.close()