from .roots.athena_analysis import *
from datetime import datetime
import gc
#from pympler.tracker import SummaryTracker
#tracker = SummaryTracker()

logging.basicConfig(filename=file.logs_loc+"/meccentricity_plot.log", encoding='utf-8', level=logging.INFO)

def main():
    aname = "_eccent" #a for analysis
    dname = "Cyl_1"
    #data_location = file.mkitp + file.cvthin2
    #data_location = file.mkitp + file.alpha3
    #data_location = "/home/mohana/Globus_data_dump/Data_backup/"
    #data_location = file.medd+"home/mohana/Globus_data_dump/Data_backup"
    #data_location = file.medd+"tmp/Data_backup"
    data_location = file.data_loc + dname
    grid_type="Cylindrical"
    savedir = file.savedir + dname + "/" + dname + aname
    mkdir_if_not_exist(savedir)
    now = datetime.now()

    for i, fnum in enumerate(range(5000, 5004, 2)):
        #2252 - 2258 scuffed
        logging.info(datetime.now()-now)
        now = datetime.now()
        logging.info("fnum = %d" % fnum)
        filename = "%s/disk.out1.%05d.athdf" % (data_location, fnum)

        aa = Athena_Analysis(filename=filename, grid_type=grid_type)

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
        vert = 2
        horz = 2
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

        #'''
        ax_mag = fig.add_subplot(gs[0, 0])
        ax_dir = fig.add_subplot(gs[0, 1])
        ax_rho = fig.add_subplot(gs[1, 1])
        ax_rhv = fig.add_subplot(gs[1, 0])

        aa.midplane_colorplot(eccent, ax_mag, log=False, vbound=[0,1], slicetype='z')
        aa.midplane_colorplot(lrl_orient, ax_dir, log=False, vbound=[0,2], angular=True, slicetype='z')
        aa.midplane_colorplot(aa.rho, ax_rho, plot_COM=True, vbound=[1e-5,1e2], slicetype='z')
        aa.midplane_colorplot(aa.rho, ax_rhv, plot_COM=False, vbound=[1e-5,1e2], slicetype='y')
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
        plt.savefig("%s/%s%s%05d.png" % (savedir, dname, aname, fnum))
        plt.close()
        del aa
        gc.collect()
        #tracker.print_diff()
        


main()
#os.system('umount ~/mnt/kitp')
os.system('umount ~/mnt/edd')