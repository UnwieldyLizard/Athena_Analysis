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
        
#main()
#os.system('umount ~/mnt/kitp')
#os.system('umount ~/mnt/edd')