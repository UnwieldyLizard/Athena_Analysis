from .roots.athena_analysis import *
import logging

logging.basicConfig(filename='/home/morgan/qresearch/profiles.log', encoding='utf-8', level=logging.INFO)

def main(starting_orbit = 35, filenum_seperation = 10, thetaslice_center = (np.pi/2), thetaslice_width = (np.pi/2), rslicepoint = 6, auto_clean=True):
    dname = "CVThick3_Profile"
    #data_location = "/home/morgan/mnt/kitp/data2/cvdisk/CVThin2/Data/"
    data_location = "/home/morgan/mnt/kitp/data/cvdisk/CVThick3/"
    savedir = "/mnt/c/Users/morga/Desktop/research_stuff/processed_data/%s" % (dname)
    mkdir_if_not_exist(savedir)

    #computing usable params based on given params
    starting_file_number = int(((starting_orbit * sim.binary_period) / sim.timesteps_per_filenum) / filenum_seperation) * filenum_seperation
    theta_start = (thetaslice_center) - (thetaslice_width/2)
    theta_end = (thetaslice_center) + (thetaslice_width/2)
    print("thetas", theta_start, theta_end)

    #initialize dictionaries
    time_avg_radial: dict[str, list] = {}
    time_avg_vert: dict[str, list] = {}
    for i, fnum in enumerate(range(starting_file_number, starting_file_number + sim.filenums_per_orbit, filenum_seperation)):
        print("fnum = %d" % fnum)
        filename = "%s/disk.out1.%05d.athdf" % (data_location, fnum)

        #debugger space
        if i == 1:
            break

        #debugger space end

        aa = Athena_Analysis(filename=filename)
        aa._axes()
        aa.get_primaries(get_rho=True, get_press=True, get_vel_r=True ,get_vel_phi=True)
        aa.get_Bfields()

        if i == 0:
            #set up and sizing
            theta_start_idx = np.argmin(abs(aa.possible_theta - theta_start))
            theta_end_idx = np.argmin(abs(aa.possible_theta - theta_end))
            r_idx = np.argmin(abs(aa.possible_r - rslicepoint))
            #print('r_idx:', r_idx, 'theta_idx:', theta_start_idx, theta_end_idx)
            angular = aa.possible_theta[theta_start_idx:theta_end_idx]
            angular = angular / np.pi
            radial = aa.possible_r
            intphitheta, intphi = aa.integrate(1, variable='phi', second_variable='theta', second_bounds=[np.pi/4, 3*np.pi/4], intermediates=True) #grid is built here 
            intphi_uncleaned = intphi.copy()
            if auto_clean == True:
                cleanup_mask = np.argwhere(intphi[theta_start_idx:theta_end_idx, r_idx] != 0).flatten()  #removes all the ghost points
                intphi = intphi[theta_start_idx:theta_end_idx, r_idx][cleanup_mask]
                angular = angular[cleanup_mask]
            else:
                intphi = intphi[theta_start_idx:theta_end_idx, r_idx]

            #building presized arrays
            time_avg_radial["surf_dens"] = np.zeros(len(aa.possible_r))
            time_avg_radial["shell_avg_press"] = np.zeros(len(aa.possible_r))
            time_avg_radial["shell_avg_Bpress"] = np.zeros(len(aa.possible_r))
            time_avg_radial["shell_avg_alpha_reyn"] = np.zeros(len(aa.possible_r))
            time_avg_radial["shell_avg_alpha_maxwell"] = np.zeros(len(aa.possible_r))
            time_avg_radial["mdot_r"] = np.zeros(len(aa.possible_r))

            time_avg_vert["azavg_dens"] = np.zeros(len(angular))
            time_avg_vert["azavg_press"] = np.zeros(len(angular))
            time_avg_vert["azavg_Bpress"] = np.zeros(len(angular))
            time_avg_vert["azavg_alpha_reyn"] = np.zeros(len(angular))
            time_avg_vert["azavg_alpha_maxwell"] = np.zeros(len(angular))
            time_avg_vert["azavg_rhovr"] = np.zeros(len(angular))

            #building constant grid attributes
            aa.get_face_areas(get_rcc_face_areas=True)    
            rccfa = aa.rcc_face_area
            grid_r = aa.r #grid was built above for first integral
            grid_theta = aa.theta
            grid_dphi = aa.dphi

        if i != 0:
            #recovering and reinserting constant grid attributes
            aa.r = grid_r
            aa.theta = grid_theta
            aa.dphi = grid_dphi

        #prep for renyold stress calculation
        shellvphi, azvphi = aa.integrate(aa.vel_phi, variable='phi', second_variable='theta', second_bounds=[np.pi/4, 3*np.pi/4], intermediates=True)
        azavgvphi = azvphi / intphi_uncleaned

        #computation of turbulent phi velocity (still part of prep for reynolds)
        turbulentvphi = np.zeros((aa.NumMeshBlocks, aa.phi_len, aa.theta_len, aa.r_len)) #defined as difference from mean phi velocity
        for n in range(aa.NumMeshBlocks):
            for i in range(aa.theta_len):
                t = np.argwhere(aa.possible_theta == aa.theta_primitive[n, i])
                for j in range(aa.r_len):
                    r = np.argwhere(aa.possible_r == aa.r_primitive[n, j])
                    turbulentvphi[n, :, i, j] = aa.vel_phi[n, :, i, j] - azavgvphi[t, r]

        #Integrals   
        shellrho, azrho = aa.integrate(aa.rho, variable='phi', second_variable='theta', second_bounds=[np.pi/4, 3*np.pi/4], intermediates=True)
        shellpress, azpress = aa.integrate(aa.press, variable='phi', second_variable='theta', second_bounds=[np.pi/4, 3*np.pi/4], intermediates=True)
        shellBpress, azBpress = aa.integrate((aa.B_r ** 2 + aa.B_theta ** 2 + aa.B_phi ** 2)/2, variable='phi', second_variable='theta', second_bounds=[np.pi/4, 3*np.pi/4], intermediates=True) 
        shellrhovr, azrhovr = aa.integrate(aa.rho * aa.vel_r, variable='phi', second_variable='theta', second_bounds=[np.pi/4, 3*np.pi/4], intermediates=True)
        shellmdot, azmdot = aa.integrate(aa.rho * aa.vel_r * rccfa, variable='phi', second_variable='theta', second_bounds=[np.pi/4, 3*np.pi/4], intermediates=True)
        shellreynolds, azreynolds = aa.integrate(aa.rho * aa.vel_r * turbulentvphi, variable='phi', second_variable='theta', second_bounds=[np.pi/4, 3*np.pi/4], intermediates=True)
        shellmaxwell, azmaxwell = aa.integrate(-aa.B_r * aa.B_phi, variable='phi', second_variable='theta', second_bounds=[np.pi/4, 3*np.pi/4], intermediates=True)

        #Shell Averages
        shellavgrho = shellrho / intphitheta
        shellavgpress = shellpress / intphitheta
        shellavgBpress = shellBpress / intphitheta
        shellavgrhovr = shellrhovr / intphitheta
        shellavgmdot = shellmdot / intphitheta
        shellavgreynolds = shellreynolds / intphitheta
        shellavgmaxwell = shellmaxwell / intphitheta

        #Azmuthal Averages
        if auto_clean == True:
            azavgrho = azrho[theta_start_idx:theta_end_idx, r_idx][cleanup_mask] / intphi
            azavgpress = azpress[theta_start_idx:theta_end_idx, r_idx][cleanup_mask] / intphi
            azavgBpress = azBpress[theta_start_idx:theta_end_idx, r_idx][cleanup_mask] / intphi
            azavgrhovr = azrhovr[theta_start_idx:theta_end_idx, r_idx][cleanup_mask] / intphi
            azavgmdot = azmdot[theta_start_idx:theta_end_idx, r_idx][cleanup_mask] / intphi
            azavgreynolds = azreynolds[theta_start_idx:theta_end_idx, r_idx][cleanup_mask] / intphi
            azavgmaxwell = azmaxwell[theta_start_idx:theta_end_idx, r_idx][cleanup_mask] / intphi
        else:
            azavgrho = azrho[theta_start_idx:theta_end_idx, r_idx] / intphi
            azavgpress = azpress[theta_start_idx:theta_end_idx, r_idx] / intphi
            azavgBpress = azBpress[theta_start_idx:theta_end_idx, r_idx] / intphi
            azavgrhovr = azrhovr[theta_start_idx:theta_end_idx, r_idx] / intphi
            azavgmdot = azmdot[theta_start_idx:theta_end_idx, r_idx] / intphi
            azavgreynolds = azreynolds[theta_start_idx:theta_end_idx, r_idx] / intphi
            azavgmaxwell = azmaxwell[theta_start_idx:theta_end_idx, r_idx] / intphi

        #computing effective alphas
        shellavg_alpha_reynolds = shellavgreynolds / shellavgpress
        azavg_alpha_reynolds = azavgreynolds / azavgpress
        shellavg_alpha_maxwell = shellavgmaxwell / shellavgpress
        azavg_alpha_maxwell = azavgmaxwell / azavgpress

        time_avg_radial["surf_dens"] += (shellavgrho / sim.binary_period)
        time_avg_radial["shell_avg_press"] += (shellavgpress / sim.binary_period)
        time_avg_radial["shell_avg_Bpress"] += (shellavgBpress / sim.binary_period)
        time_avg_radial["shell_avg_alpha_reyn"] += (shellavg_alpha_reynolds / sim.binary_period)
        time_avg_radial["shell_avg_alpha_maxwell"] += (shellavg_alpha_maxwell / sim.binary_period)
        time_avg_radial["mdot_r"] += (shellavgmdot / sim.binary_period)

        time_avg_vert["azavg_dens"] += (azavgrho / sim.binary_period)
        time_avg_vert["azavg_press"] += (azavgpress / sim.binary_period)
        time_avg_vert["azavg_Bpress"] += (azavgBpress / sim.binary_period)
        time_avg_vert["azavg_alpha_reyn"] += (azavg_alpha_reynolds / sim.binary_period)
        time_avg_vert["azavg_alpha_maxwell"] += (azavg_alpha_maxwell / sim.binary_period)
        time_avg_vert["azavg_rhovr"] += (azavgrhovr / sim.binary_period)

    print("\tplot")
    # radial profile
    vert = 4
    horz = 1
    gs = gridspec.GridSpec(vert, horz)
    fig = plt.figure(figsize=(horz*3, vert*2), dpi=300)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(radial, time_avg_radial["surf_dens"], "C0-")
    ax.set_yscale("log")
    ax.set_title("surface density (log scale)")
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    ax = fig.add_subplot(gs[1, 0])
    ax.plot(radial, time_avg_radial["shell_avg_press"], "C0-", label="gas")
    ax.plot(radial, time_avg_radial["shell_avg_Bpress"], "C1-", label="B")
    ax.set_yscale("log")
    ax.set_title("pressure (log scale)")
    ax.legend(loc="upper right")
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    ax = fig.add_subplot(gs[2, 0])
    ax.plot(radial, time_avg_radial["shell_avg_alpha_reyn"], "C0-", label="reyn")
    ax.plot(radial, time_avg_radial["shell_avg_alpha_maxwell"], "C1-", label="B")
    ax.set_yscale("log")
    ax.set_ylim([1e-4, 1])
    ax.set_title(r"$\alpha$ (log scale)")
    ax.legend(loc="upper right")
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    ax = fig.add_subplot(gs[3, 0])
    ax.plot(radial, time_avg_radial["mdot_r"], "C0-")
    ax.set_xlabel("r")
    #ax.set_yscale("log")
    ax.set_title(r"$\dot M$")
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    plt.tight_layout()
    plt.savefig("%s/%s_radial.png" % (savedir, dname))
    plt.close()

    # vertical profile
    vert = 4
    horz = 1
    gs = gridspec.GridSpec(vert, horz)
    fig = plt.figure(figsize=(horz*3, vert*2), dpi=300)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(angular, time_avg_vert["azavg_dens"], "C0-")
    #ax.set_ylabel("density")
    ax.set_yscale("log")
    ax.set_title("volume density (log scale)")
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    ax = fig.add_subplot(gs[1, 0])
    ax.plot(angular, time_avg_vert["azavg_press"], "C0-", label="gas")
    ax.plot(angular, time_avg_vert["azavg_Bpress"], "C1-", label="B")
    #ax.set_ylabel("density")
    ax.set_yscale("log")
    ax.set_title("pressure (log scale)")
    ax.legend(loc="upper right")
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    ax = fig.add_subplot(gs[2, 0])
    ax.plot(angular, time_avg_vert["azavg_alpha_reyn"], "C0-", label="reyn")
    ax.plot(angular, time_avg_vert["azavg_alpha_maxwell"], "C1-", label="B")
    #ax.set_ylabel("density")
    ax.set_yscale("log")
    ax.set_ylim([1e-4, 1])
    ax.set_title(r"$\alpha$ (log scale)")
    ax.legend(loc="upper right")
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    ax = fig.add_subplot(gs[3, 0])
    ax.plot(angular, time_avg_vert["azavg_rhovr"], "C0-")
    ax.set_xlabel(r"$\theta$")
    #ax.set_ylabel("density")
    ax.set_title(r"$\rho v_r$")
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    plt.tight_layout()
    plt.savefig("%s/%s_angular.png" % (savedir, dname))
    plt.close()
    
    '''
    # save all data
    with open("%s/%s_pickle.dat" % (savedir, dname), "wb") as pickle_file:
        pickle.dump({
            "angular": angular,
            "time_avg_vert": time_avg_vert,
            "radial": radial,
            "time_avg_radial": time_avg_radial,
        }, pickle_file)
    '''   

def profile_movie(coordinates, filenum_seperation = 2, vertslice_center = (np.pi/2), rslicepoint = 6, auto_clean=True, MAGNETIC_FIELDS_ENABLED=True):
    aname = "_profile" #a for analysis
    #dname = "CVThick3_Profile"
    dname = "Cyl_1"
    #data_location = "/home/morgan/mnt/kitp/data2/cvdisk/CVThin2/Data/"
    #data_location = "/home/morgan/mnt/kitp/data/cvdisk/CVThick3/"
    data_location = file.data_loc + "Data_backup"
    savedir = file.savedir + dname + "/" + dname + aname
    mkdir_if_not_exist(savedir)

    if MAGNETIC_FIELDS_ENABLED:
        final_terms = ["rho", "press", "Bpress", "alpha_reyn", "alpha_maxwell"]
        terms = ["rho", "press", "rhovr", "mdot", "reynolds", "Bpress", "maxwell"]
    else:
        final_terms = ["rho", "press", "alpha_reyn"]
        terms = ["rho", "press", "rhovr", "mdot", "reynolds"]
    
    #computing usable params based on given params
    if coordinates == "Spherical":
        thetaslice_width = (np.pi/2)
        theta_start = (vertslice_center) - (thetaslice_width/2)
        theta_end = (vertslice_center) + (thetaslice_width/2)
    
    loops_per_orbit = int(sim.filenums_per_orbit / filenum_seperation)
    now = datetime.now()

    for i, fnum in enumerate(range(0, 10000, filenum_seperation)):
        logging.info(datetime.now()-now)
        now = datetime.now()
        logging.info("fnum = %d" % fnum)
        filename = "%s/disk.out1.%05d.athdf" % (data_location, fnum)    

        aa = Athena_Analysis(filename=filename, grid_type=coordinates)
        aa._axes()
        aa.get_primaries(get_rho=True, get_press=True, get_vel_r=True ,get_vel_phi=True)
        if MAGNETIC_FIELDS_ENABLED:
            aa.get_Bfields()

        if i == 0:
            #set up and sizing
            aa.native_grid()
            r_idx = np.argmin(abs(aa.possible_r - rslicepoint))
            radial = aa.possible_r
            if coordinates == "Spherical":
                theta_start_idx = np.argmin(abs(aa.possible_theta - theta_start))
                theta_end_idx = np.argmin(abs(aa.possible_theta - theta_end))
                vert = aa.possible_theta[theta_start_idx:theta_end_idx]
                vert = vert / np.pi
                intshell, intphi = aa.integrate(1, variable="shell", intermediates=True)
                intphi_uncleaned = intphi.copy()
                if auto_clean == True:
                    cleanup_mask = np.argwhere(intphi[theta_start_idx:theta_end_idx, r_idx] != 0).flatten()  #removes all the ghost points
                    intphi = intphi[theta_start_idx:theta_end_idx, r_idx][cleanup_mask]
                    vert = vert[cleanup_mask]
                else:
                    intphi = intphi[theta_start_idx:theta_end_idx, r_idx]
            elif coordinates == "Cylindrical":
                vert = aa.possible_z
                intshell, intphi = aa.integrate(1, variable='shell', intermediates=True) #grid is built here
                intphi_uncleaned = intphi.copy()
                if auto_clean == True:
                    cleanup_mask = np.argwhere(intphi[:, r_idx] != 0).flatten()  #removes all the ghost points
                    intphi = intphi[:, r_idx][cleanup_mask]
                    vert = vert[cleanup_mask]
                else:
                    intphi = intphi[:, r_idx] 
            r_idx = np.argmin(abs(aa.possible_r - rslicepoint))
            radial = aa.possible_r

            #initialize dictionaries
            time_avg_radial: dict[str, list] = {}
            time_avg_vert: dict[str, list] = {}
            radials = np.array([{}] * loops_per_orbit)
            verts = np.array([{}] * loops_per_orbit)

            #old terms below for reference
            '''
            time_avg_radial["surf_dens"] = np.zeros(len(aa.possible_r))
            time_avg_radial["shell_avg_press"] = np.zeros(len(aa.possible_r))
            time_avg_radial["shell_avg_alpha_reyn"] = np.zeros(len(aa.possible_r))
            time_avg_radial["mdot_r"] = np.zeros(len(aa.possible_r))
            time_avg_vert["azavg_dens"] = np.zeros(len(vert))
            time_avg_vert["azavg_press"] = np.zeros(len(vert))
            time_avg_vert["azavg_alpha_reyn"] = np.zeros(len(vert))
            time_avg_vert["azavg_rhovr"] = np.zeros(len(vert))
            if MAGNETIC_FIELDS_ENABLED:
                time_avg_radial["shell_avg_Bpress"] = np.zeros(len(aa.possible_r))
                time_avg_radial["shell_avg_alpha_maxwell"] = np.zeros(len(aa.possible_r))
                time_avg_vert["azavg_Bpress"] = np.zeros(len(vert))
                time_avg_vert["azavg_alpha_maxwell"] = np.zeros(len(vert))
            '''
            for j in range(loops_per_orbit):
                radials[j] = {}
                verts[j] = {}
                for key in final_terms:
                    radials[i]["shell_avg_"+key] = np.zeros(len(radial))
                    verts[i]["azavg_"+key] = np.zeros(len(vert))
                radials[i]["mdot_r"] = np.zeros(len(radial))
                verts[i]["azavg_rhovr"] = np.zeros(len(vert))

            #building constant grid attributes
            aa.get_face_areas(get_rcc_face_areas=True)    
            rccfa = aa.rcc_face_area
            grid_r = aa.r #grid was built above for first integral
            if coordinates == "Spherical":
                grid_theta = aa.theta
            if coordinates == "Cylindrical":
                grid_z = aa.z
            grid_dphi = aa.dphi

        if i != 0:
            #recovering and reinserting constant grid attributes
            aa.r = grid_r
            if coordinates == "Spherical":
                aa.theta = grid_theta
            if coordinates == "Cylindrical":
                aa.z = grid_z
            aa.dphi = grid_dphi

        #prep for renyold stress calculation
        shellvphi, azvphi = aa.integrate(aa.vel_phi, variable="shell", intermediates=True)
        azavgvphi = azvphi / intphi_uncleaned

        #computation of turbulent phi velocity (still part of prep for reynolds)
        turbulentvphi = np.zeros(aa.array_size) #defined as difference from mean phi velocity
        if coordinates == "Spherical":
            for n in range(aa.NumMeshBlocks):
                for k in range(aa.theta_len):
                    t = np.argwhere(aa.possible_theta == aa.theta_primitive[n, k])
                    for j in range(aa.r_len):
                        r = np.argwhere(aa.possible_r == aa.r_primitive[n, j])
                        turbulentvphi[n, :, k, j] = aa.vel_phi[n, :, k, j] - azavgvphi[t, r]
        if coordinates == "Cylindrical":
            for n in range(aa.NumMeshBlocks):
                for k in range(aa.z_len):
                    z = np.argwhere(aa.possible_z == aa.z_primitive[n, k])
                    for j in range(aa.r_len):
                        r = np.argwhere(aa.possible_r == aa.r_primitive[n, j])
                        turbulentvphi[n, k, :, j] = aa.vel_phi[n, k, :, j] - azavgvphi[z, r]

        #Integrals
        shell = {} #for shell integrals
        az = {} #for azmuthal integrals
        shell["rho"], az["rho"] = aa.integrate(aa.rho, variable='shell', intermediates=True)
        shell["press"], az["press"] = aa.integrate(aa.press, variable='shell', intermediates=True)
        shell["rhovr"], az["rhovr"] = aa.integrate(aa.rho * aa.vel_r, variable='shell', intermediates=True)
        shell["mdot"], az["mdot"] = aa.integrate(aa.rho * aa.vel_r * rccfa, variable='shell', intermediates=True)
        shell["reynolds"], az["reynolds"] = aa.integrate(aa.rho * aa.vel_r * turbulentvphi, variable='shell', intermediates=True)
        if MAGNETIC_FIELDS_ENABLED:
            shell["Bpress"], az["Bpress"] = aa.integrate((aa.B_r ** 2 + aa.B_theta ** 2 + aa.B_phi ** 2)/2, variable='shell', intermediates=True)
            shell["maxwell"], az["maxwell"] = aa.integrate(-aa.B_r * aa.B_phi, variable='shell', intermediates=True)

        #Shell Averages
        shellavg = {} # for shell averages
        for key in shell:
            shellavg[key] = shell[key] / intshell

        #Azmuthal Averages
        azavg = {} # for azmuthal averages
        if coordinates == "Spherical":
            if auto_clean == True:
                for key in az:
                    azavg[key] = az[key][theta_start_idx:theta_end_idx, r_idx][cleanup_mask] / intphi
            else:
                for key in az:
                    azavg[key] = az[key][theta_start_idx:theta_end_idx, r_idx] / intphi
        elif coordinates == "Cylindrical":
            if auto_clean == True:
                for key in az:
                    azavg[key] = az[key][:, r_idx][cleanup_mask] / intphi
            else:
                for key in az:
                    azavg[key] = az[key][:, r_idx] / intphi

        #computing effective alphas
        shellavg["alpha_reyn"] = shellavg["reynolds"] / shellavg["press"]
        azavg["alpha_reyn"] = azavg["reynolds"] / azavg["press"]
        if MAGNETIC_FIELDS_ENABLED:
            shellavg["alpha_maxwell"]= shellavg["maxwell"] / shellavg["press"]
            azavg["alpha_maxwell"] = azavg["maxwell"] / azavg["press"]

        if i < loops_per_orbit:
            for key in final_terms:
                radials[i]["shell_avg_"+key] = shellavg[key]
                verts[i]["azavg_"+key] = azavg[key]
            radials[i]["mdot_r"] = shellavg["mdot"]
            verts[i]["azavg_rhovr"] = azavg["rhovr"]
        else:
            radials[:-1] = radials[1:]
            verts[:-1] = verts[1:]
            for key in final_terms:
                radials[-1]["shell_avg_"+key] = shellavg[key]
                verts[-1]["azavg_"+key] = azavg[key]
            radials[-1]["mdot_r"] = shellavg["mdot"]
            verts[-1]["azavg_rhovr"] = azavg["rhovr"]

        #sizing arrays and wiping old data
        for key in final_terms:
            time_avg_radial["shell_avg_"+key] = np.zeros(len(radial))
            time_avg_vert["azavg_"+key] = np.zeros(len(vert))
        time_avg_radial["mdot_r"] = np.zeros(len(radial))
        time_avg_vert["azavg_rhovr"] = np.zeros(len(vert))

        if i < loops_per_orbit:
            for l in range(i+1): #orginally divided by sim.binary_period
                for key in radials[0]:
                    time_avg_radial[key] += radials[l][key] / (i+1)
                for key in verts[0]:
                    time_avg_vert[key] += verts[l][key] / (i+1)
        else:
            for l in range(loops_per_orbit): #orginally divided by sim.binary_period
                for key in radials[0]:
                    time_avg_radial[key] += radials[l][key] / loops_per_orbit
                for key in verts[0]:
                    time_avg_vert[key] += verts[l][key] / loops_per_orbit

        logging.info("\tplot")
        # radial profile
        grid_vert = 4
        grid_horz = 2
        gs = gridspec.GridSpec(grid_vert, grid_horz)
        fig = plt.figure(figsize=(grid_horz*3, grid_vert*2), dpi=300)
        fig.suptitle(f'Orbit: {(aa.time/sim.binary_period):.2f}')

        ax = fig.add_subplot(gs[0, 0])
        ax.plot(radial, time_avg_radial["shell_avg_rho"], "C0-")
        ax.set_yscale("log")
        ax.set_ylim([1e-4, 5e1])
        ax.set_title("surface density (log scale)")
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
        #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

        ax = fig.add_subplot(gs[1, 0])
        ax.plot(radial, time_avg_radial["shell_avg_press"], "C0-", label="gas")
        if MAGNETIC_FIELDS_ENABLED:
            ax.plot(radial, time_avg_radial["shell_avg_Bpress"], "C1-", label="B")
        ax.set_yscale("log")
        ax.set_ylim([1e-4, 1e3])
        ax.set_title("pressure (log scale)")
        ax.legend(loc="upper right")
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
        #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

        ax = fig.add_subplot(gs[2, 0])
        ax.plot(radial, time_avg_radial["shell_avg_alpha_reyn"], "C0-", label="reyn")
        if MAGNETIC_FIELDS_ENABLED:
            ax.plot(radial, time_avg_radial["shell_avg_alpha_maxwell"], "C1-", label="B")
        ax.set_yscale("log")
        ax.set_ylim([1e-4, 1e2])
        ax.set_title(r"$\alpha$ (log scale)")
        ax.legend(loc="upper right")
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
        #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

        ax = fig.add_subplot(gs[3, 0])
        ax.plot(radial, time_avg_radial["mdot_r"], "C0-")
        ax.set_xlabel("r")
        #ax.set_yscale("log")
        ax.set_ylim([-0.005, 0.005])
        ax.set_title(r"$\dot M$")
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
        #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

        '''
        plt.tight_layout()
        plt.savefig("%s/%s_profile_radial%05d.png" % (savedir, dname, fnum))
        plt.close()

        # vertical profile
        grid_vert = 4
        grid_horz = 1
        gs = gridspec.GridSpec(grid_vert, grid_horz)
        fig = plt.figure(figsize=(grid_horz*3, grid_vert*2), dpi=300)
        fig.suptitle(f'Orbit: {(aa.time/sim.binary_period):.2f}')
        '''

        ax = fig.add_subplot(gs[0, 1])
        ax.plot(vert, time_avg_vert["azavg_rho"], "C0-")
        #ax.set_ylabel("density")
        ax.set_yscale("log")
        ax.set_ylim([1e-4, 5e1])
        ax.set_title("volume density (log scale)")
        if coordinates == "Spherical":
            ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
        else:
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

        ax = fig.add_subplot(gs[1, 1])
        ax.plot(vert, time_avg_vert["azavg_press"], "C0-", label="gas")
        if MAGNETIC_FIELDS_ENABLED:
            ax.plot(vert, time_avg_vert["azavg_Bpress"], "C1-", label="B")
        #ax.set_ylabel("density")
        ax.set_yscale("log")
        ax.set_ylim([1e-4, 1e3])
        ax.set_title("pressure (log scale)")
        ax.legend(loc="upper right")
        if coordinates == "Spherical":
            ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
        else:
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

        ax = fig.add_subplot(gs[2, 1])
        ax.plot(vert, time_avg_vert["azavg_alpha_reyn"], "C0-", label="reyn")
        if MAGNETIC_FIELDS_ENABLED:
            ax.plot(vert, time_avg_vert["azavg_alpha_maxwell"], "C1-", label="B")
        #ax.set_ylabel("density")
        ax.set_yscale("log")
        ax.set_ylim([1e-4, 1e2])
        ax.set_title(r"$\alpha$ (log scale)")
        ax.legend(loc="upper right")
        if coordinates == "Spherical":
            ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
        else:
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

        ax = fig.add_subplot(gs[3, 1])
        ax.plot(vert, time_avg_vert["azavg_rhovr"], "C0-")
        #ax.set_ylabel("density")
        ax.set_ylim([-0.75, 0.25])
        ax.set_title(r"$\rho v_r$")
        if coordinates == "Spherical":
            ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
            ax.set_xlabel(r"$\theta$")
        else:
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
            ax.set_xlabel("z")

        '''
        plt.tight_layout()
        plt.savefig("%s/%s_profile_vert%05d.png" % (savedir, dname, fnum))
        plt.close()
        '''
        plt.tight_layout()
        plt.savefig("%s/%s%s%05d.png" % (savedir, dname, aname, fnum))
        plt.close()

#main(auto_clean=True)
profile_movie(coordinates="Cylindrical", vertslice_center=0, MAGNETIC_FIELDS_ENABLED=False)