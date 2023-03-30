#Eventually needs to contain rewrite of mprofile code
from .roots.athena_analysis import *

logging.basicConfig(filename=file.logs_loc+"/profile.log", encoding='utf-8', level=logging.INFO)

def mass_profile(dname, fnum, r_slicepoint, grid_type):
    aname = "_profile" #a for analysis
    data_location = file.data_loc + dname
    grid_type=grid_type
    savedir = file.savedir + dname + "/" + dname + aname
    mkdir_if_not_exist(savedir)
    filename = "%s/disk.out1.%05d.athdf" % (data_location, fnum)
    aa = Athena_Analysis(filename=filename, grid_type=grid_type)

    aa.get_primaries(get_rho=True)
    aa.native_grid()
    r_idx = np.argmin(abs(aa.possible_r - r_slicepoint))
    rad_rho, az_rho = aa.integrate(aa.rho, "shell", intermediates=True)
    rad_rho = rad_rho / aa.integrate(1, "shell")
    vert_rho = az_rho[:,r_idx] / (2*np.pi*aa.possible_r[r_idx])

    vert = 1
    horz = 2
    gs = gridspec.GridSpec(vert, horz)
    fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

    vertical_distribution = np.arange(0,40,1)
    
    ax_rad = fig.add_subplot(gs[0, 0])
    ax_vert = fig.add_subplot(gs[0, 1])

    ax_rad.plot(aa.possible_r, rad_rho)
    ax_rad.plot(np.full(len(vertical_distribution), sim.three_one_res), vertical_distribution, "C3--", label="Resonant Radius")
    #ax_rad.plot(np.full(len(vertical_distribution), sim.old_three_one_res), vertical_distribution, "C1--", label="Old Resonant Radius")
    if grid_type=="Cylindrical":
        ax_vert.plot(aa.possible_z, vert_rho)
    if grid_type=="Spherical":
        ax_vert.plot(aa.possible_theta, vert_rho)

    ax_rad.set_ylim([0,40])
    ax_vert.set_ylim([0,40])
    ax_rad.set_title(r"Radial $\rho$")
    ax_vert.set_title(r"Vertical $\rho$"+" (r=%s)" % (r_slicepoint))
    ax_rad.legend()

    plt.tight_layout()
    orbit = (aa.time / sim.binary_period)
    fig.suptitle(f"orbit: {orbit:.2f}")
    plt.savefig("%s/%s_density%s%05d.png" % (savedir, dname, aname, fnum))
    plt.close()

def mass_profile_loop(dname, fnum_range, file_spacing, r_slicepoint, grid_type):
    now = datetime.now()
    fnum_range = np.arange(fnum_range[0], fnum_range[-1]+1, file_spacing)
    for fnum in fnum_range:
        logging.info(datetime.now()-now)
        now = datetime.now()
        logging.info("fnum = %d" % fnum)
        mass_profile(dname, fnum, r_slicepoint, grid_type)