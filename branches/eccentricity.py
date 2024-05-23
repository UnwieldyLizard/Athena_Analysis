from .roots.athena_analysis import *
from .roots.misc_func import *
from datetime import datetime
import gc
#from pympler.tracker import SummaryTracker
#tracker = SummaryTracker()

logging.basicConfig(filename=file.logs_loc+"/eccentricity.log", encoding='utf-8', level=logging.INFO)

class Eccentricity:
    def __init__(self, dname, start_point = 0):
        self.dname = dname
        self.aname = "_eccent"
        self.file_spacing = find_file_spacing(dname, start_point)
        self.data_location = file.data_loc + dname
        self.grid_type = file.grid_types[dname]
        self.savedir = file.savedir + dname + "/" + dname + self.aname
        mkdir_if_not_exist(self.savedir)

    def plot_loop(self, fnum_range, inertial=True):
        now = datetime.now()

        for i, fnum in enumerate(range(fnum_range[0], fnum_range[1], file_spacing)):
            #2252 - 2258 scuffed
            logging.info(datetime.now()-now)
            now = datetime.now()
            logging.info("fnum = %d" % fnum)
            filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)

            aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)

            if i == 0:
                aa.native_grid(get_r=True, get_phi=True)
                r_grid = aa.r
                phi_grid = aa.phi
                phi_len = aa.phi_len
                if grid_type == "Spherical":
                    theta_len = aa.theta_len
                if grid_type == "Cylindrical":
                    z_len = aa.z_len
                r_len = aa.r_len
                aa._build_angular_cmap()
                ang_cmap = aa.angular_cmap
                vec_ar_siz = aa.vector_array_size
                ar_siz = aa.array_size

                """
                #initializing averages/sums
                if grid_type == "Spherical":
                    sum_lrl_orient = np.zeros((aa.NumMeshBlocks, aa.phi_len, aa.theta_len, aa.r_len))
                    sum_eccent = np.zeros((aa.NumMeshBlocks, aa.phi_len, aa.theta_len, aa.r_len))
                    ave_lrl_orient = np.zeros((aa.NumMeshBlocks, aa.phi_len, aa.theta_len, aa.r_len))
                    ave_eccent = np.zeros((aa.NumMeshBlocks, aa.phi_len, aa.theta_len, aa.r_len))
                """

            if i != 0:
                aa.r = r_grid
                aa.phi = phi_grid
                aa.phi_len = phi_len
                if grid_type == "Spherical":
                    aa.theta_len = theta_len
                if grid_type == "Cylindrical":
                    aa.z_len = z_len
                aa.r_len = r_len
                aa.angular_cmap = ang_cmap
                aa.vector_array_size = vec_ar_siz
                aa.array_size = ar_siz

            aa.get_primaries(get_rho=True, get_vel_r=True, get_vel_phi=True)

            arctan_start = datetime.now()
            eccent, lrl_orient = aa.get_lrl()
            '''
            #computing sums
            sum_lrl_orient += lrl_orient
            sum_eccent += eccent

            #computing averages
            ave_lrl_orient = (sum_lrl_orient / (i+1))
            ave_eccent = (sum_eccent / (i+1))
            '''

            #Transforming to inertial frame
            if inertial:
                rotation = -1*sim.orbital_Omega * aa.time
                lrl_orient = (lrl_orient - rotation) % (2*np.pi)
            else:
                rotation = 0

            vert = 2
            horz = 2
            gs = gridspec.GridSpec(vert, horz)
            fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

            #'''
            ax_mag = fig.add_subplot(gs[0, 0])
            ax_dir = fig.add_subplot(gs[0, 1])
            ax_rho = fig.add_subplot(gs[1, 1])
            ax_rhv = fig.add_subplot(gs[1, 0])

            aa.midplane_colorplot(eccent, ax_mag, log=False, vbound=[0,1], slicetype='z', rotation=rotation)
            aa.midplane_colorplot(lrl_orient, ax_dir, log=False, vbound=[0,2], angular=True, slicetype='z', rotation=rotation)
            aa.midplane_colorplot(aa.rho, ax_rho, plot_COM=False, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
            aa.midplane_colorplot(aa.rho, ax_rhv, plot_COM=False, vbound=[1e-5,1e2], slicetype='vert', rotation=rotation)
            ax_mag.set_title("Eccentricity Magnitude")
            ax_dir.set_title("LRL Angle from X Axis")
            ax_rho.set_title("Density Midplane Slice")
            ax_rhv.set_title("Density Profile Slice")
            '''
            ax_x = fig.add_subplot(gs[0, 0])
            ax_y = fig.add_subplot(gs[0, 1])
            ax_r = fig.add_subplot(gs[1, 0])
            ax_phi = fig.add_subplot(gs[1, 1])

            aa.midplane_colorplot(lrl_x, ax_x, log=False, vbound=[-2,2])
            aa.midplane_colorplot(lrl_y, ax_y, log=False, vbound=[-2,2])
            aa.midplane_colorplot(lrl_r, ax_r, log=False, vbound=[-2,2])
            aa.midplane_colorplot(lrl_phi, ax_phi, log=False, vbound=[-2,2])
            #aa.midplane_colorplot(lrl_orient, ax_y, log=False, vbound=[0,2], angular=True)
            ax_x.set_title("LRL x component")
            ax_y.set_title("LRL y component")
            ax_r.set_title("LRL r component")
            ax_phi.set_title("LRL phi component")
            #ax_y.set_title("LRL Angle from X Axis")
            '''

            plt.tight_layout()
            orbit = (aa.time / sim.binary_period)
            plt.subplots_adjust(top=(1-0.01*(16/vert)))
            fig.suptitle(f"orbit: {orbit:.2f}")
            plt.savefig("%s/%s%s%05d.png" % (self.savedir, self.dname, self.aname, fnum))
            plt.close()
            del aa
            gc.collect()
            #tracker.print_diff()
    
    def plot(self, fnum, inertial=True, sname=""):
        logging.info("fnum = %d" % fnum)
        filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)

        aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)

        aa.get_primaries(get_rho=True, get_vel_r=True, get_vel_phi=True)

        eccent, lrl_orient = aa.get_lrl()
        '''
        #computing sums
        sum_lrl_orient += lrl_orient
        sum_eccent += eccent

        #computing averages
        ave_lrl_orient = (sum_lrl_orient / (i+1))
        ave_eccent = (sum_eccent / (i+1))
        '''

        #Transforming to inertial frame
        if inertial:
            rotation = -1*sim.orbital_Omega * aa.time
            lrl_orient = (lrl_orient - rotation) % (2*np.pi)
        else:
            rotation = 0

        vert = 2
        horz = 2
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

        #'''
        ax_mag = fig.add_subplot(gs[0, 0])
        ax_dir = fig.add_subplot(gs[0, 1])
        ax_rho = fig.add_subplot(gs[1, 1])
        ax_rhv = fig.add_subplot(gs[1, 0])

        aa.midplane_colorplot(eccent, ax_mag, log=True, vbound=[1e-2,1], slicetype='z', rotation=rotation, sci_notation=False)     
        aa.midplane_colorplot(lrl_orient, ax_dir, log=False, vbound=[0,2], angular=True, slicetype='z', rotation=rotation)
        aa.midplane_colorplot(aa.rho, ax_rho, plot_COM=False, vbound=[1e-5,1e2], slicetype='z', rotation=rotation)
        aa.midplane_colorplot(aa.rho, ax_rhv, plot_COM=False, vbound=[1e-5,1e2], slicetype='vert', rotation=rotation)
        ax_mag.set_title("Eccentricity Magnitude")
        ax_dir.set_title("LRL Angle from X Axis")
        ax_rho.set_title("Density Midplane Slice")
        ax_rhv.set_title("Density Profile Slice")
        '''
        ax_x = fig.add_subplot(gs[0, 0])
        ax_y = fig.add_subplot(gs[0, 1])
        ax_r = fig.add_subplot(gs[1, 0])
        ax_phi = fig.add_subplot(gs[1, 1])

        aa.midplane_colorplot(lrl_x, ax_x, log=False, vbound=[-2,2])
        aa.midplane_colorplot(lrl_y, ax_y, log=False, vbound=[-2,2])
        aa.midplane_colorplot(lrl_r, ax_r, log=False, vbound=[-2,2])
        aa.midplane_colorplot(lrl_phi, ax_phi, log=False, vbound=[-2,2])
        #aa.midplane_colorplot(lrl_orient, ax_y, log=False, vbound=[0,2], angular=True)
        ax_x.set_title("LRL x component")
        ax_y.set_title("LRL y component")
        ax_r.set_title("LRL r component")
        ax_phi.set_title("LRL phi component")
        #ax_y.set_title("LRL Angle from X Axis")
        '''

        plt.tight_layout()
        orbit = (aa.time / sim.binary_period)
        plt.subplots_adjust(top=(1-0.01*(16/vert)))
        fig.suptitle(f"{file.display_name[self.dname]}  Orbit: {orbit:.2f}")
        plt.savefig("%s/%s%s%05d%s.png" % (self.savedir, self.dname, self.aname, fnum, sname))
        plt.close()
        del aa
        gc.collect()
        #tracker.print_diff()

    def two_point_comparison(self, start_fnum=None, start_time=None, end_fnum=None, delta_t=None):
        if start_fnum is None and start_time is None:
            raise("You must proved either an explicit end point or a start time")
        if start_fnum is None:
            start_fnum = int(sim.filenums_per_orbit * start_time)
        if end_fnum is None and delta_t is None:
            raise("You must proved either an explicit end point or a duration")
        if end_fnum is None:
            end_fnum = int(start_fnum + sim.filenums_per_orbit * delta_t)

        sname = f"Comp{start_fnum}-{end_fnum}"

        start_point = Athena_Analysis(filename="%s/disk.out1.%05d.athdf" % (self.data_location, start_fnum), grid_type=self.grid_type)
        end_point = Athena_Analysis(filename="%s/disk.out1.%05d.athdf" % (self.data_location, end_fnum), grid_type=self.grid_type)

        start_point.get_primaries(get_rho=True, get_vel_r=True, get_vel_phi=True)
        end_point.get_primaries(get_rho=True, get_vel_r=True, get_vel_phi=True)

        start_eccent, start_lrl_orient = start_point.get_lrl()
        end_eccent, end_lrl_orient = end_point.get_lrl()

        diff = end_eccent - start_eccent

        int_diff = start_point.integrate(diff * start_point.rho * end_point.rho, "shell") / start_point.integrate(start_point.rho * end_point.rho, "shell")
  
        vert = 1
        horz = 1
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

        #'''
        ax_mag = fig.add_subplot(gs[0, 0])
        
        ax_mag.plot(start_point.possible_r, int_diff)
        #start_point.midplane_colorplot(end_eccent-start_eccent, ax_mag, log=True, vbound=[-1,1], slicetype='z', sci_notation=False)     
        ax_mag.set_title(r"$\Delta$e between orbits"+f" {start_point.time/sim.binary_period:.3} - {end_point.time/sim.binary_period:.3}")
        
        plt.tight_layout()
        plt.subplots_adjust(top=(1-0.01*(16/vert)))
        fig.suptitle(f"{file.display_name[self.dname]}")
        plt.savefig("%s/%s%s%s.png" % (self.savedir, self.dname, self.aname, sname))
        plt.close()
        gc.collect()
        #tracker.print_diff()

    def flux(self, start_fnum, fnum_len=1):
        end_fnum = start_fnum + fnum_len - 1

        for i, fnum in enumerate(range(start_fnum, end_fnum+1, 1)):
            filename = "%s/disk.out1.%05d.athdf" % (self.data_location, fnum)

            aa = Athena_Analysis(filename=filename, grid_type=self.grid_type)
            
            aa.get_primaries(get_vel_r=True, get_rho=True)

            if i == 0:
                aa.native_grid(get_r=True)
                r_axis = aa.possible_r
                flux = np.zeros(len(aa.possible_r))
                M_dot = np.zeros(len(aa.possible_r))

            lrl_native, lrl_cart = aa.get_lrl(components = True)    

            flux += aa.integrate(aa.vel_r * aa.rho * vec.get_magnitude(lrl_cart), "Shell") / (fnum_len)
            M_dot += aa.integrate(aa.vel_r * aa.rho, "Shell") /(fnum_len)

        # growth plot
        vert = 2
        horz = 1
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(2*horz*3, vert*3), dpi=300)
        
        ax = fig.add_subplot(gs[0, 0])
        #ax.plot(r_axis, flux)
        ax.plot(r_axis, flux)
        #ax.legend()
        ax.set_xlabel("r")
        ax.set_ylabel("Flux")
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

        ax = fig.add_subplot(gs[1, 0])
        #ax.plot(r_axis, flux)
        ax.plot(r_axis, M_dot)
        #ax.legend()
        ax.set_xlabel("r")
        ax.set_ylabel(r"$\dot{M}$")
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

        plt.subplots_adjust(top=(1-0.01*(16/vert)))
        start_orbit = start_fnum / sim.filenums_per_orbit
        end_orbit = end_fnum / sim.filenums_per_orbit
        if start_fnum == end_fnum:
            title = file.display_name[self.dname]+f": orbit {start_orbit:.2f}"
        else:
            title = file.display_name[self.dname]+f": orbits {start_orbit:.2f} - {end_orbit:.2f}" 
        plt.suptitle(title)
        plt.tight_layout()
        plt.savefig("%s/%s_flux_%05d-%05d.png" % (self.savedir, self.dname, start_fnum, end_fnum))
        plt.close()

        
#main()
#os.system('umount ~/mnt/kitp')
#os.system('umount ~/mnt/edd')