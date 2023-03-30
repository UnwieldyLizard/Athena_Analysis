from .roots.athena_analysis import *

logging.basicConfig(filename=file.logs_loc+"/orbital.log", encoding='utf-8', level=logging.INFO)

def orbital_velocity_analysis(dname, fnum, grid_type, MHD=False, alpha=None, phi_slicepoint=None):
    aname = "_orbital" #a for analysis
    data_location = file.data_loc + dname
    grid_type=grid_type
    savedir = file.savedir + dname + "/" + dname + aname
    mkdir_if_not_exist(savedir)
    filename = "%s/disk.out1.%05d.athdf" % (data_location, fnum)
    aa = Athena_Analysis(filename=filename, grid_type=grid_type)
    
    aa.native_grid(get_r=True, get_z=True, get_phi=True)
    aa.get_primaries(get_rho=True, get_press=True, get_vel_r=True, get_vel_phi=True, get_vel_z=True)
    aa.get_potentials(get_accel=True, get_companion_grav=True)
    if phi_slicepoint is None:
        phi_idx = "az_ave"
    else:
        phi_idx = np.argmin(abs(aa.possible_phi - phi_slicepoint))
        [phi_slice_x, phi_slice_y] = [aa.possible_r * np.cos(aa.possible_phi[phi_idx]), aa.possible_r * np.sin(aa.possible_phi[phi_idx])]        

    v_phi = reduce_dimensions(aa.vel_phi, phi_idx, aa=aa)
    v_r = reduce_dimensions(aa.vel_r, phi_idx, aa=aa)
    rad_rho = reduce_dimensions(aa.rho, phi_idx, aa=aa)
    #rad_press = aa.integrate(aa.press, "shell") / shell_normalization

    accel_list = ["White Dwarf Gravity", "Companion Gravity", "Radial Flow", "Pressure"]
    accels = {}
    accels["White Dwarf Gravity"] = -1*(sim.gm1 / (aa.possible_r ** 2))
    accels["Companion Gravity"] = reduce_dimensions(-1*aa.differentiate(aa.accel_pot + aa.companion_grav_pot, "r"), phi_idx, aa=aa)
    accels["Pressure"] = reduce_dimensions(-1 *aa.differentiate(aa.press, "r") / aa.rho, phi_idx, aa=aa)
    if MHD:
        aa.get_Bfields()
        accel_list.append("Magnetic Pressure")
        accel_list.append("Magnetic Tension")
        accels["Magnetic Pressure"] = (reduce_dimensions(
            -1* aa.differentiate((aa.B_z**2+aa.B_phi**2+aa.B_r**2), "r") / aa.rho, phi_idx, aa=aa))
        accels["Magnetic Tension"] = (reduce_dimensions(
            ((-1*aa.B_phi*aa.B_phi/aa.r)
            + aa.B_z * aa.differentiate(aa.B_r, 'z')
            + (aa.B_phi/aa.r) * aa.differentiate(aa.B_r, 'phi')
            + aa.B_r * aa.differentiate(aa.B_r, 'r'))
            / aa.rho, phi_idx, aa=aa))
    else:
        accels_list.append("Viscosity")
        accels["Viscosity"] = reduce_dimensions(((-3/2) * alpha * (aa.differentiate(aa.press, "phi")/aa.r)) / aa.rho, phi_idx, aa=aa)
    """
    #companion gravity, the long hard way
    beta = (aa.r / sim.bin_sep)
    laplace_expression = 1 / (1-2*beta*np.cos(aa.phi)+(beta**2))**(1/2)
    laplace_coefficient = (1/np.pi)*aa.integrate(laplace_expression, variable='z', second_variable='phi') / aa.possible_r #normalized in z cux height is 1
    d_laplace_coefficient = np.zeros(len(aa.possible_r))
    d_laplace_coefficient[0] = 0
    d_laplace_coefficient[1:] = (laplace_coefficient[1:] - laplace_coefficient[:-1])
    d_beta = (aa.possible_dr_primitive / sim.bin_sep)
    laplace_factor = d_laplace_coefficient / d_beta
    accels["Companion Gravity"] = aa.possible_r*(((sim.gm2)/(2*(sim.bin_sep**2)*aa.possible_r))*(laplace_factor))
    """
    #radial acceleration
    i = 1
    found_previous_file = False
    while found_previous_file == False:
        try:
            filename_last = filename = "%s/disk.out1.%05d.athdf" % (data_location, fnum - i)
            aa_last = Athena_Analysis(filename=filename_last, grid_type=grid_type)
            found_previous_file = True
        except:
            i += 1
    print(fnum-i)
    aa_last.get_primaries(get_vel_r=True)
    v_r_last = reduce_dimensions(aa_last.vel_r, phi_idx, aa)
    dv_rdt = (v_r - v_r_last) / (i*sim.timesteps_per_filenum)
    matderiv_r = reduce_dimensions((#(-1*aa.vel_phi*aa.vel_phi/aa.r) 
        + aa.vel_z * aa.differentiate(aa.vel_r, 'z')
        #+ aa.vel_phi * (aa.differentiate(aa.vel_r, 'phi')/aa.r)
        + aa.vel_r * aa.differentiate(aa.vel_r, 'r')), phi_idx, aa=aa)
    a_r = dv_rdt + matderiv_r #note this is not the full acceleration field, just the components which do not depend on v_phi. See theory pdf
    accels["Radial Flow"] = -1*a_r
    
    mat_deriv_phi_coef = reduce_dimensions((aa.differentiate(aa.vel_r, 'phi')/aa.r), phi_idx, aa=aa)
    total_accel = np.zeros(len(aa.possible_r))
    for key in accel_list:
        total_accel += accels[key]
    predicted_v_phi = ((aa.possible_r/2)*(
        + mat_deriv_phi_coef
        - 2*sim.orbital_Omega
        + np.sqrt((mat_deriv_phi_coef**2) - 4*sim.orbital_Omega*mat_deriv_phi_coef - (4/aa.possible_r)*(total_accel))))

    #plotting
    vert = 2
    horz = 4
    gs = gridspec.GridSpec(vert, horz)
    fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

    ax_rho = fig.add_subplot(gs[0, 1])
    ax_rad = fig.add_subplot(gs[0, 0])
    ax_vphi = fig.add_subplot(gs[1, 0])
    ax_vphi_log = fig.add_subplot(gs[1, 1])
    ax_accel = fig.add_subplot(gs[0, 2])
    ax_accel_log = fig.add_subplot(gs[0, 3])
    ax_accels = fig.add_subplot(gs[1, 2])
    ax_accels_log = fig.add_subplot(gs[1, 3])
    
    #mins = np.array([min(predicted_v_phi), min(v_phi)])
    #maxs = np.array([max(predicted_v_phi), max(v_phi)])
    vertical_distribution = np.arange(min([min(v_phi),0]), max(v_phi), 1)

    ax_rad.plot(aa.possible_r, rad_rho, "b", label=r"$\rho$")
    #ax_rad.plot(aa.possible_r, rad_press, c="orange", label="P")
    ax_rad.plot(np.full(len(vertical_distribution), sim.three_one_res), vertical_distribution, "C3--", label="Resonant Radius")
    ax_rad.set_ylim([0,40])
    ax_rad.legend()
    aa.midplane_colorplot(aa.rho, ax_rho, plot_COM=False, vbound=[1e-5, 1e2], slicetype='z')
    if phi_idx != "az_ave":
        ax_rho.plot(phi_slice_x, phi_slice_y, "r")
    ax_rho.set_title(r"$\rho$")

    ax_vphi.plot(aa.possible_r, v_phi, "C2-", label=r"$v_{\phi}$")
    ax_vphi.plot(aa.possible_r, predicted_v_phi, "b--", label=r"Predicted $v_{\phi}$")
    ax_vphi.plot(np.full(len(vertical_distribution), sim.three_one_res), vertical_distribution, "C3--", label="Res Radius")
    ax_vphi.plot(aa.possible_r, 2*sim.orbital_Omega*aa.possible_r, "C1--", label="Res Velocity")
    ax_vphi.set_ylim([min(v_phi), max(v_phi)])
    ax_vphi_log.plot(aa.possible_r, v_phi, "C2-", label=r"$v_{\phi}$")
    ax_vphi_log.plot(aa.possible_r, predicted_v_phi, "b--",label=r"Predicted $v_{\phi}$")
    ax_vphi_log.plot(np.full(len(vertical_distribution), sim.three_one_res), vertical_distribution, "C3--", label="Res Radius")
    ax_vphi_log.plot(aa.possible_r, 2*sim.orbital_Omega*aa.possible_r, "C1--", label="Res Velocity")
    ax_vphi_log.set_xscale("log")
    ax_vphi_log.set_yscale("log")
    ax_vphi.legend()
    ax_vphi_log.legend()

    for k, key in enumerate(accel_list):
        ax_accels.plot(aa.possible_r, -1*accels[key], f"C{k}-", label=key)
        ax_accels_log.plot(aa.possible_r, -1*accels[key], f"C{k}-", label=key)
        ax_accels_log.set_xscale("log")
        ax_accels_log.set_yscale("log")
    ax_accels.legend()
    ax_accels.set_ylim([-300,300])
    ax_accels_log.legend()
    
    ax_accel.plot(aa.possible_r, ((v_phi**2)/aa.possible_r), "C2-", label="Measured")
    ax_accel.plot(aa.possible_r, -1*total_accel, "C9--", label="Total Source Terms")
    ax_accel_log.plot(aa.possible_r, ((v_phi**2)/aa.possible_r), "C2-", label="Measured")
    ax_accel_log.plot(aa.possible_r, -1*total_accel, "C9--", label="Total Source Terms")
    ax_accel_log.set_xscale("log")
    ax_accel_log.set_yscale("log")
    ax_accel.set_title("Centripetal Acceleration")
    ax_accel_log.set_title("Centripetal Acceleration (Log)")
    ax_accel.legend()
    ax_accel_log.legend()

    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    orbit = (aa.time / sim.binary_period)
    if phi_idx == "az_ave":
        fig.suptitle(f"orbit: {orbit:.2f}, azimuthal ave")
        plt.savefig("%s/%s%s%05d_az_ave.png" % (savedir, dname, aname, fnum))
    else:
        fig.suptitle(f"orbit: {orbit:.2f}, phi={aa.possible_phi[phi_idx]/np.pi:.2f}pi")
        plt.savefig("%s/%s%s%05d_phi=" % (savedir, dname, aname, fnum) +f"{aa.possible_phi[phi_idx]/np.pi:.2f}pi"+".png")
    plt.close()

def reduce_dimensions(q, phi_idx, aa):
    if phi_idx == "az_ave":
        if not hasattr(reduce_dimensions, "normalization_weight"):
            reduce_dimensions.normalization_weight = aa.integrate(1, "shell")
            print("this should print exactly once")
        return (aa.integrate(q, "shell") / reduce_dimensions.normalization_weight)
    else:
        if (not hasattr(reduce_dimensions, "normalization_weight")) and aa.gridtype == "Spherical":
            reduce_dimensions.normalization_weight = aa.integrate(1, "theta")
            print("this should print exactly once")
        if aa.gridtype == "Spherical":
            return((aa.integrate(q, "theta") / reduce_dimensions.normalization_weight)[phi_idx])
        if (not hasattr(reduce_dimensions, "normalization_weight")) and aa.gridtype == "Cylindrical":
            reduce_dimensions.normalization_weight = aa.integrate(1, "z")
            print("this should print exactly once")
        if aa.gridtype == "Cylindrical":
            return((aa.integrate(q, "z") / reduce_dimensions.normalization_weight)[phi_idx])

def orbital_velocity_analysis_loop(dname, fnum_range, file_spacing, grid_type):
    now = datetime.now()
    fnum_range = np.arange(fnum_range[0], fnum_range[-1]+1, file_spacing)
    for fnum in fnum_range:
        logging.info(datetime.now()-now)
        now = datetime.now()
        logging.info("fnum = %d" % fnum)
        orbital_velocity_analysis(dname, fnum, grid_type)

def res_orbit_plot(dname, fnum, grid_type, ax=None):
    aname = "_resonant_orbit" #a for analysis
    data_location = file.data_loc + dname
    grid_type=grid_type
    savedir = file.savedir + dname + "/" + dname + aname
    mkdir_if_not_exist(savedir)
    filename = "%s/disk.out1.%05d.athdf" % (data_location, fnum)
    aa = Athena_Analysis(filename=filename, grid_type=grid_type)
    
    aa.native_grid(get_r=True, get_z=True, get_phi=True)
    aa.get_primaries(get_rho=True, get_vel_phi=True)

    res_radii = np.zeros(len(aa.possible_phi))
    v_phi = aa.integrate(aa.vel_phi, "z")
    res_vel = aa.integrate(2*sim.orbital_Omega*aa.r, "z")
    for p in range(len(aa.possible_phi)):
        r_idx = np.argmin(abs(v_phi[p] - res_vel[p]))
        res_radii[p] = aa.possible_r[r_idx]
    
    if ax is None:
        vert = 1
        horz = 1
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
        ax = fig.add_subplot(gs[0, 0])
        auto_plot = True

    aa.midplane_colorplot(aa.rho, ax=ax)
    ax.plot(res_radii * np.cos(aa.possible_phi), res_radii * np.sin(aa.possible_phi), "b" ,label="Resonant Radius", linewidth=0.5)
    #ax.legend(loc='upper right', bbox_to_anchor=(1.5, 0.5))

    if auto_plot:
        plt.savefig("%s/%s%s%05d.png" % (savedir, dname, aname, fnum))
        plt.close()