from .roots.athena_analysis import *
import logging
import pickle

logging.basicConfig(filename=file.logs_loc+"/meccentricity_profile.log", encoding='utf-8', level=logging.INFO)


def main(dname, fnum_limits, file_spacing, plot_every=100, pickle_every=None, MHD = True, coordinates = 'Spherical', pickle_flags=[], restart=False, alpha=None):
    aname = "_eccent_growth_prec" #a for analysis
    data_location = file.data_loc + dname
    #data_location = file.mkitp + file.alpha3
    savedir = file.savedir + dname + "/" + dname + aname
    mkdir_if_not_exist(savedir)
    mkdir_if_not_exist(savedir + "/pickles")
    if pickle_every == None:
        pickle_every = plot_every
    fnum_range = list(range(fnum_limits[0], fnum_limits[1], file_spacing))
    pickle_flags = np.array(pickle_flags)
    pickle_flags = np.append(pickle_flags, fnum_limits[1])
    if restart == True:
        dname = dname + "_%ssrt" % fnum_range[0]

    if MHD == True:
        terms = ["tidal", "press", "boundary", "Bpress", "Btens", "vpress", "hpress", "vBpress", "hBpress"]
        sum_terms = ["tidal", "press", "boundary", "Bpress", "Btens"]
        MAGNETIC_FIELDS_ENABLED = True
    else:
        terms = ["tidal", "press", "boundary", "visc", "vpress", "hpress", "tidal_vr"]
        sum_terms = ["tidal", "press", "boundary", "visc"]
        MAGNETIC_FIELDS_ENABLED = False
    
    orbit_series = np.zeros(len(fnum_range))
    time_series = np.zeros(len(fnum_range))
    eccent_series = np.zeros(len(fnum_range))
    eccent_phase_series = np.zeros(len(fnum_range))

    eccent_term_series = {}
    integrated_eccent_term_series = {}
    eccent_term_phase_series = {}
    integrated_eccent_term_phase_series = {}

    for key in terms:
        eccent_term_series[key] = np.zeros(len(fnum_range))
        integrated_eccent_term_series[key] = np.zeros(len(fnum_range))
        integrated_eccent_term_series[key][0] = 0
        eccent_term_phase_series[key] = np.zeros(len(fnum_range))
        integrated_eccent_term_phase_series[key] = np.zeros(len(fnum_range))
        integrated_eccent_term_phase_series[key][0] = 0
    integrated_phase_flips = np.zeros(len(fnum_range))

    eccent_sum_term_series = np.zeros(len(fnum_range))
    integrated_eccent_sum_term_series = np.zeros(len(fnum_range))
    integrated_eccent_sum_term_series[0] = 0
    eccent_phase_sum_term_series = np.zeros(len(fnum_range))
    integrated_eccent_phase_sum_term_series = np.zeros(len(fnum_range))
    integrated_eccent_phase_sum_term_series[0] = 0

    now = datetime.now()
    for i, fnum in enumerate(fnum_range):

        if (i == 0) and (fnum != 0) and (restart == False):
            logging.info("setting up")
            found_load_point = False
            load_point = fnum - file_spacing
            while found_load_point == False:
                if os.path.exists("%s/pickles/%s_pickle_%05d.dat" % (savedir, dname, load_point)):
                    logging.info("Found data, loading from: %s" % load_point)
                    found_load_point = True
                else:
                    load_point -= file_spacing
                if load_point <= 0:
                    raise("No load point, you need to restart")
            with open("%s/pickles/%s_pickle_%05d.dat" % (savedir, dname, load_point), "rb") as pickle_file:
                data = pickle.load(pickle_file)
            i_0 = data["offset"]

            orbit_series = np.zeros(len(fnum_range)+i_0)
            time_series = np.zeros(len(fnum_range)+i_0)
            eccent_series = np.zeros(len(fnum_range)+i_0)
            eccent_phase_series = np.zeros(len(fnum_range)+i_0)

            logging.info(str(i_0))
            logging.info(str(len(fnum_range)))

            eccent_term_series = {}
            integrated_eccent_term_series = {}
            eccent_term_phase_series = {}
            integrated_eccent_term_phase_series = {}

            for key in terms:
                eccent_term_series[key] = np.zeros(len(fnum_range)+i_0)
                integrated_eccent_term_series[key] = np.zeros(len(fnum_range)+i_0)
                integrated_eccent_term_series[key][0] = 0
                eccent_term_phase_series[key] = np.zeros(len(fnum_range)+i_0)
                integrated_eccent_term_phase_series[key] = np.zeros(len(fnum_range)+i_0)
                integrated_eccent_term_phase_series[key][0] = 0
            integrated_phase_flips = np.zeros(len(fnum_range)+i_0)

            eccent_sum_term_series = np.zeros(len(fnum_range)+i_0)
            integrated_eccent_sum_term_series = np.zeros(len(fnum_range)+i_0)
            integrated_eccent_sum_term_series[0] = 0
            eccent_phase_sum_term_series = np.zeros(len(fnum_range)+i_0)
            integrated_eccent_phase_sum_term_series = np.zeros(len(fnum_range)+i_0)
            integrated_eccent_phase_sum_term_series[0] = 0                

            time_series[:i_0+1] = data["time_series"]
            orbit_series[:i_0+1] = data["orbit_series"]
            eccent_series[:i_0+1] = data["eccent_series"]
            eccent_phase_series[:i_0+1] = data["eccent_phase_series"]
            for key in terms:
                eccent_term_series[key][:i_0+1] = data["eccent_term_series"][key][:i_0+1]
                eccent_term_phase_series[key][:i_0+1] = data["eccent_term_phase_series"][key][:i_0+1]
                integrated_eccent_term_series[key][:i_0+1] = data["integrated_eccent_term_series"][key][:i_0+1]
                integrated_eccent_term_phase_series[key][:i_0+1] = data["integrated_eccent_term_phase_series"][key][:i_0+1]
            eccent_sum_term_series[:i_0+1] = data["eccent_sum_term_series"]
            eccent_phase_sum_term_series[:i_0+1] = data["eccent_phase_sum_term_series"]
            integrated_eccent_sum_term_series[:i_0+1] = data["integrated_eccent_sum_term_series"]
            integrated_eccent_phase_sum_term_series[:i_0+1] = data["integrated_eccent_phase_sum_term_series"]
            integrated_phase_flips[:i_0+1] = data["integrated_phase_flips"][:i_0+1]

            old_mass_weighted_eccent = data["old_mass_weighted_eccent"]
            
            #logging.info("tidal_int: "+ " ".join(map(str, [integrated_eccent_term_series["tidal"][i_0+1]])) )

            logging.info("set up done")
        elif i == 0:
            i_0 = -1

        j = i + i_0 + 1     

        logging.info(datetime.now()-now)
        now = datetime.now()

        logging.info("fnum = %d" % fnum)
        logging.info("i: "+ " ".join(map(str, [i])) + " j: " + " ".join(map(str, [j]))) 
        filename = "%s/disk.out1.%05d.athdf" % (data_location, fnum)
        
        
        aa = Athena_Analysis(filename=filename, grid_type= coordinates)
        
        orbit_series[j] = (aa.time / sim.binary_period)
        time_series[j] = (aa.time)

        #calculate forces
        aa.get_primaries(get_rho=True, get_press=True)
        aa.get_potentials(get_companion_grav=True, get_accel=True, get_wd_grav=True)
        aa.build_vectors()
        if MAGNETIC_FIELDS_ENABLED:
            aa.get_Bfields()
        force = {}
        force["tidal"] = -1 * aa.rho * aa.gradient(aa.accel_pot + aa.companion_grav_pot, coordinates=coordinates)
        force["press"] = -1 * aa.gradient(aa.press, coordinates=coordinates)
        if MAGNETIC_FIELDS_ENABLED:
            if aa.gridtype == "Spherical":
                force["Bpress"] = -1 * aa.gradient(((aa.B_r ** 2) + (aa.B_theta ** 2) + (aa.B_phi ** 2)) / 2, coordinates=coordinates)
                force["Btens"] = aa.material_derivative([aa.B_phi, aa.B_theta, aa.B_r], [aa.B_phi, aa.B_theta, aa.B_r])
                force["vBpress"] = np.array(force["Bpress"], dtype=np.float32)
                force["vBpress"][0, :] *= 0
                force["vBpress"][2, :] *= 0
            elif aa.gridtype == "Cylindrical":
                force["Bpress"] = -1 * aa.gradient(((aa.B_z ** 2) + (aa.B_phi ** 2) + (aa.B_r ** 2)) / 2, coordinates=coordinates)
                force["Btens"] = aa.material_derivative([aa.B_z, aa.B_phi, aa.B_r], [aa.B_z, aa.B_phi, aa.B_r])
                force["vBpress"] = np.array(force["Bpress"], dtype=np.float32)
                force["vBpress"][1, :] *= 0
                force["vBpress"][2, :] *= 0
            force["hBpress"] = force["Bpress"] - force["vBpress"]
        else:
            force["visc"] = aa.alpha_visc(alpha)
            logging.info("alpha visc in testing")
        force["vpress"] = np.array(force["press"], dtype=np.float32)
        if aa.gridtype == "Spherical":
            force["vpress"][0, :] *= 0
            force["vpress"][2, :] *= 0
        elif aa.gridtype == "Cylindrical":
            force["vpress"][1, :] *= 0
            force["vpress"][2, :] *= 0
        force["hpress"] = force["press"] - force["vpress"]

        #Mass Weighted Eccentricity
        lrl_native, lrl_cart = aa.get_lrl(components = True)        
        total_mass = aa.integrate(aa.rho, variable='all')
        rhoxlrl = np.array([aa.rho*lrl_cart[0], aa.rho*lrl_cart[1], np.zeros(aa.array_size)])
        total_rhoxlrl = aa.vec_vol_integrate(rhoxlrl)
        mass_weighted_eccent = total_rhoxlrl / total_mass

        test_cart = aa.vec_vol_integrate(lrl_cart)

        '''
        logging.info("integrateed lrl cart: \n" + " ".join(map(str, test_cart)))
        logging.info("mass weighted eccent: \n" + " ".join(map(str, mass_weighted_eccent)))
        logging.info("mass weighted eccent magnitude: \n" + " ".join(map(str, [vec.get_magnitude(mass_weighted_eccent)])))
        '''

        eccent_series[j] = vec.get_magnitude(mass_weighted_eccent)
        eccent_phase_series[j] = np.arctan2(mass_weighted_eccent[1], mass_weighted_eccent[0]) + sim.orbital_Omega * aa.time

        # calculate C
        aa.get_primaries(get_vel_r=True, get_vel_phi=True)
        C = {}
        for key in force:
            C[key] = aa.get_C_vec(force[key])
            C[key] = aa.native_to_cart(C[key])
            evol_C = aa.vec_vol_integrate(C[key]) / total_mass
            growth, prec = vec.growth_prec(evol_C, mass_weighted_eccent)

            eccent_term_series[key][j] = growth
            eccent_term_phase_series[key][j] = prec

            #logging.info("%s: \n" % key + " ".join(map(str, evol_C)))

        #Boundary
        rhoxlrl_flux = np.array([aa.get_boundary_flux(rhoxlrl[0]), aa.get_boundary_flux(rhoxlrl[1]), 0])
        advect_evol = -1 * rhoxlrl_flux / total_mass
        adv_growth, adv_prec = vec.growth_prec(advect_evol, mass_weighted_eccent)

        mass_flux = aa.get_boundary_flux(aa.rho)

        massd_evol = mass_flux * np.array(mass_weighted_eccent) / total_mass
        growth, prec = vec.growth_prec(massd_evol, mass_weighted_eccent)

        eccent_term_series["boundary"][j] = adv_growth + growth
        eccent_term_phase_series["boundary"][j] = adv_prec + prec


        '''
        logging.info("J (should be same):" + " ".join(map(str, [j])))
        logging.info("this growth rate: "+" ".join(map(str, [eccent_term_series["tidal"][j]])))
        logging.info("last growth rate: "+" ".join(map(str, [eccent_term_series["tidal"][j-1]])))
        logging.info("ave growth rate: "+" ".join(map(str, [(eccent_term_series["tidal"][j] + eccent_term_series["tidal"][j-1])/2])))
        logging.info("time delta: "+" ".join(map(str, [(time_series[j]-time_series[j-1])])))
        logging.info("eccent growth: "+" ".join(map(str, [(time_series[j]-time_series[j-1])*(eccent_term_series["tidal"][j] + eccent_term_series["tidal"][j-1])/2])))
        '''

        #integrate values forward
        total_eccent_increment = 0
        total_phase_increment = 0
        for key in terms:
            if j != 0:
                eccent_increment = (time_series[j]-time_series[j-1])*(eccent_term_series[key][j] + eccent_term_series[key][j-1])/2
                integrated_eccent_term_series[key][j] = integrated_eccent_term_series[key][j-1] + eccent_increment
                phase_increment = (time_series[j]-time_series[j-1])*(eccent_term_phase_series[key][j] + eccent_term_phase_series[key][j-1])/2
                integrated_eccent_term_phase_series[key][j] = integrated_eccent_term_phase_series[key][j-1] + phase_increment

                total_eccent_increment += eccent_increment
                total_phase_increment += phase_increment

        #summing all contributions
        for key in sum_terms:
            eccent_sum_term_series[j] += eccent_term_series[key][j]
            eccent_phase_sum_term_series[j] += eccent_term_phase_series[key][j]
        
        if j == 0 or j == 1:
            integrated_eccent_sum_term_series[i] = eccent_series[i]
            integrated_eccent_phase_sum_term_series[i] = eccent_phase_series[i]
        elif j != 0 and j != 1:
            integrated_eccent_sum_term_series[j] = integrated_eccent_sum_term_series[j-1] + total_eccent_increment
            #checking for flip and adjusting
            if integrated_eccent_sum_term_series[j] < 0:
                flip_phase = np.pi
                integrated_eccent_sum_term_series[j] = np.abs(integrated_eccent_sum_term_series[j])
                logging.info("auto-flipped phase")
            elif (vec.ortho_dot(old_mass_weighted_eccent, mass_weighted_eccent)) < 0:
                flip_phase = np.pi
                logging.info("manually flipped phase")
            else:
                flip_phase = 0
            integrated_phase_flips[j] = integrated_phase_flips[j-1] + flip_phase
            integrated_eccent_phase_sum_term_series[j] = integrated_eccent_phase_sum_term_series[j-1] + total_phase_increment + flip_phase
        
        #memorising old eccentricity vector for comparison in next time step
        old_mass_weighted_eccent = mass_weighted_eccent

        '''
        #manual readout
        logging.info("sum eccent: \n" + " ".join(map(str, [integrated_eccent_sum_term_series[j]])))
        logging.info("sum prec: \n" + " ".join(map(str, [integrated_eccent_phase_sum_term_series[j]])))
        for key in terms:
            #logging.info("%s eccent: \n" % key + " ".join(map(str, eccent_term_series[key][j])))
            #logging.info("%s prec: \n" % key + " ".join(map(str, eccent_term_phase_series[key][j])))
            logging.info("%s integrated eccent: \n" % key + " ".join(map(str, [integrated_eccent_term_series[key][j]])))
            logging.info("%s integrated prec: \n" % key + " ".join(map(str, [integrated_eccent_term_phase_series[key][j]])))
        '''
        
        # set up plot
        if j % plot_every == 0 or i == 0:
            logging.info("\tplot")

            # growth plot
            vert = 2
            horz = 1
            gs = gridspec.GridSpec(vert, horz)
            fig = plt.figure(figsize=(2*horz*3, vert*3), dpi=300)

            ax = fig.add_subplot(gs[0, 0])
            ax.plot(orbit_series[:j+1], eccent_series[:j+1], "C2-", label="measured")
            ax.plot(orbit_series[:j+1], integrated_eccent_sum_term_series[:j+1], "C9--", label="total source terms")
            ax.set_xlabel("binary orbit")
            ax.set_ylabel("eccent magnitude")
            ax.set_ylim([1e-4, 1])
            ax.set_yscale("log")
            #ax.set_title("mhd_3d eccent magnitude")
            ax.legend(loc="upper left")
            #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

            ax = fig.add_subplot(gs[1, 0])
            for k, key in enumerate(terms):
                ax.plot(orbit_series[:j+1], integrated_eccent_term_series[key][:j+1], f"C{k}-", label=key)
            ax.set_xlabel("binary orbit")
            ax.set_ylabel("eccent contribution")
            #ax.set_title("mhd_3d time integrated sources")
            ax.legend(loc="upper left")
            #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

            plt.tight_layout()
            plt.savefig("%s/%s_growth_%05d.png" % (savedir, dname, fnum))
            plt.close()

            # precession plot
            vert = 2
            horz = 1
            gs = gridspec.GridSpec(vert, horz)
            fig = plt.figure(figsize=(2*horz*3, vert*3), dpi=300)

            ax = fig.add_subplot(gs[0, 0])
            ax.plot(orbit_series[:j+1], wrap_phase(eccent_phase_series[:j+1]) / np.pi, "C3-", label="measured")
            ax.plot(orbit_series[:j+1], (wrap_phase(integrated_eccent_phase_sum_term_series[:j+1])) / np.pi, "C1--", label="total source terms")
            ax.set_xlabel("binary orbit")
            ax.set_ylabel("eccent phase")
            ax.set_ylim([-1, 1])
            #ax.set_yscale("log")
            ax.set_title("mhd_3d eccent phase")
            ax.legend(loc="upper left")
            #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
            #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

            ax = fig.add_subplot(gs[1, 0])
            for k, key in enumerate(terms):
                ax.plot(orbit_series[:j+1], (wrap_phase(integrated_eccent_term_phase_series[key][:j+1])) / np.pi, f"C{k}-", label=key)
            ax.plot(orbit_series[:j+1], (wrap_phase(integrated_phase_flips[:j+1])) / np.pi, f"C3--", label="phase inversion")
            ax.set_xlabel("binary orbit")
            ax.set_ylabel("eccent phase contribution")
            ax.set_ylim([-1, 1])
            ax.set_title("mhd_3d time integrated sources")
            ax.legend(loc="upper left")
            #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
            #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

            plt.tight_layout()
            plt.savefig("%s/%s_prec_%05d.png" % (savedir, dname, fnum))
            plt.close()

        if (j % pickle_every == 0) or (fnum in pickle_flags):
            logging.info("pickling up to %s" % fnum)
            with open("%s/pickles/%s_pickle_%05d.dat" % (savedir, dname, fnum), "wb") as pickle_file:
                pickle.dump({
                    "offset": j,
                    "time_series": time_series[:j+1],
                    "orbit_series": orbit_series[:j+1],
                    "eccent_series": eccent_series[:j+1],
                    "eccent_phase_series": eccent_phase_series[:j+1],
                    "eccent_term_series": eccent_term_series,
                    "eccent_term_phase_series": eccent_term_phase_series,
                    "eccent_sum_term_series": eccent_sum_term_series[:j+1],
                    "eccent_phase_sum_term_series": eccent_phase_sum_term_series[:j+1],
                    "integrated_eccent_term_series": integrated_eccent_term_series,
                    "integrated_eccent_term_phase_series": integrated_eccent_term_phase_series,
                    "integrated_eccent_sum_term_series": integrated_eccent_sum_term_series[:j+1],
                    "integrated_eccent_phase_sum_term_series": integrated_eccent_phase_sum_term_series[:j+1],
                    "integrated_phase_flips": integrated_phase_flips[:j+1],
                    "old_mass_weighted_eccent": old_mass_weighted_eccent,
                }, pickle_file)


def replot(dname, fnum, pname="", aspect_ratio=2, MHD=False):
    """
    replots data

    Parameters
    ----------
    dname, str:
        The name of the data set
    fnum, int:
        The filenumber where you'd like to replot
    pname, str:
        Special prefix for output file names
    aspect_ratio, float:
        Length over height
    MHD, bool:
        If it's MHD or not
    """
    aname = "_eccent_growth_prec" #a for analysis
    savedir = file.savedir + dname + "/" + dname + aname
    logging.info("looking for file")
    found_load_point = False
    load_point = fnum
    while found_load_point == False:
        if os.path.exists("%s/pickles/%s_pickle_%05d.dat" % (savedir, dname, load_point)):
            logging.info("Found data, plotting from: %s" % load_point)
            found_load_point = True
        else:
            load_point -= 1
        if load_point <= 0:
            raise("No data found to plot")
    with open("%s/pickles/%s_pickle_%05d.dat" % (savedir, dname, load_point), "rb") as pickle_file:
        data = pickle.load(pickle_file)

    if MHD == True:
        terms = ["tidal", "press", "boundary", "Bpress", "Btens", "vpress", "hpress", "vBpress", "hBpress"]
        sum_terms = ["tidal", "press", "boundary", "Bpress", "Btens"]
        MAGNETIC_FIELDS_ENABLED = True
    else:
        terms = ["tidal", "press", "boundary", "visc", "vpress", "hpress", "tidal_vr"]
        sum_terms = ["tidal", "press", "boundary", "visc"]
        MAGNETIC_FIELDS_ENABLED = False
    
    # growth plot
    vert = 2
    horz = 1

    if aspect_ratio >= 1:
        vert_scale = 1
        horz_scale = aspect_ratio
    else:
        vert_scale = 1 / aspect_ratio
        horz_scale = 1

    logging.info("Aspect ratio: "+str([horz_scale, vert_scale]))

    gs = gridspec.GridSpec(vert, horz)
    fig = plt.figure(figsize=(horz_scale * horz*3, vert_scale * vert*3), dpi=300)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(data["orbit_series"], data["eccent_series"], "C2-", label="measured")
    ax.plot(data["orbit_series"], data["integrated_eccent_sum_term_series"], "C9--", label="total source terms")
    ax.set_xlabel("binary orbit")
    ax.set_ylabel("eccent magnitude")
    ax.set_ylim([1e-4, 1])
    ax.set_yscale("log")
    #ax.set_title("mhd_3d eccent magnitude")
    ax.legend(loc="upper left")
    #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    ax = fig.add_subplot(gs[1, 0])
    for k, key in enumerate(terms):
        ax.plot(data["orbit_series"], data["integrated_eccent_term_series"][key][:data["offset"]+1], f"C{k}-", label=key)
    ax.set_xlabel("binary orbit")
    ax.set_ylabel("eccent contribution")
    #ax.set_title("mhd_3d time integrated sources")
    ax.legend(loc="upper left")
    #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    plt.tight_layout()
    plt.savefig("%s/%s%s_growth_%05d.png" % (savedir, pname, dname, fnum))
    plt.close()

    # precession plot
    gs = gridspec.GridSpec(vert, horz)
    fig = plt.figure(figsize=(horz_scale*horz*3, vert_scale * vert*3), dpi=300)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(data["orbit_series"], wrap_phase(data["eccent_phase_series"]) / np.pi, "C3-", label="measured")
    ax.plot(data["orbit_series"], (wrap_phase(data["integrated_eccent_phase_sum_term_series"])) / np.pi, "C1--", label="total source terms")
    ax.set_xlabel("binary orbit")
    ax.set_ylabel("eccent phase")
    ax.set_ylim([-1, 1])
    #ax.set_yscale("log")
    ax.set_title("mhd_3d eccent phase")
    ax.legend(loc="upper left")
    #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    ax = fig.add_subplot(gs[1, 0])
    for k, key in enumerate(terms):
        ax.plot(data["orbit_series"], (wrap_phase(data["integrated_eccent_term_phase_series"][key][:data["offset"]+1])) / np.pi, f"C{k}-", label=key)
    ax.plot(data["orbit_series"], (wrap_phase(data["integrated_phase_flips"])) / np.pi, f"C3--", label="phase inversion")
    ax.set_xlabel("binary orbit")
    ax.set_ylabel("eccent phase contribution")
    ax.set_ylim([-1, 1])
    ax.set_title("mhd_3d time integrated sources")
    ax.legend(loc="upper left")
    #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    plt.tight_layout()
    plt.savefig("%s/%s%s_prec_%05d.png" % (savedir, pname, dname, fnum))
    plt.close()

def tidal_profile(dname, fnum, grid_type, phi_slicepoint=None, radial_slice_loop=False):
    aname = "_tidal" #a for analysis
    tidal_profile.aname = aname
    data_location = file.data_loc + dname
    grid_type=grid_type
    savedir = file.savedir + dname + "/" + dname + aname
    mkdir_if_not_exist(savedir)
    filename = "%s/disk.out1.%05d.athdf" % (data_location, fnum)
    aa = Athena_Analysis(filename=filename, grid_type=grid_type)
    
    aa.native_grid(get_r=True, get_z=True, get_phi=True)
    aa.get_primaries(get_rho=True, get_vel_phi=True)

    if phi_slicepoint is None:
        phi_idx = "az_ave"
    else:
        phi_idx = np.argmin(abs(aa.possible_phi - phi_slicepoint))
        [phi_slice_x, phi_slice_y] = [aa.possible_r * np.cos(aa.possible_phi[phi_idx]), aa.possible_r * np.sin(aa.possible_phi[phi_idx])]

    aa.get_potentials(get_accel=True, get_companion_grav=True)
    tidal_C = reduce_dimensions(aa.native_to_cart(aa.get_C_vec(-1 * aa.rho * aa.gradient(aa.accel_pot + aa.companion_grav_pot, coordinates=grid_type))), phi_idx, aa)
    lrl, lrl_cart = aa.get_lrl(components=True)
    eccent = reduce_dimensions(aa.rho * lrl_cart, phi_idx, aa)
    tidal_growth = (vec.ortho_dot(eccent, tidal_C) / vec.get_magnitude(eccent)) / aa.integrate(aa.rho, "all")

    vert = 1
    horz = 2
    gs = gridspec.GridSpec(vert, horz)
    fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
    ax_rho = fig.add_subplot(gs[0, 1])
    ax_tid = fig.add_subplot(gs[0, 0])
    
    ax_tid.plot(aa.possible_r, tidal_growth, "b")
    ax_tid.set_title("Tidal Eccentricity Contribution")
    ax_tid.set_xlabel("r")
    if radial_slice_loop:
        ax_tid.set_ylim([-1e-4, 1e-4])
    else:
        ax_tid.set_ylim([-1e-4, 1e-4])
    aa.midplane_colorplot(aa.rho, ax_rho, plot_COM=False, vbound=[1e-5, 1e2], slicetype='z')
    if phi_idx != "az_ave":
        ax_rho.plot(phi_slice_x, phi_slice_y, "r")
    ax_rho.set_title(r"$\rho$")

    plt.tight_layout()
    plt.subplots_adjust(top=0.82)
    orbit = (aa.time / sim.binary_period)
    if phi_idx == "az_ave":
        fig.suptitle(f"orbit: {orbit:.2f}, azimuthal ave")
        plt.savefig("%s/%s%s%05d_az_ave.png" % (savedir, dname, aname, fnum))
    else:
        fig.suptitle(f"orbit: {orbit:.2f}, phi={aa.possible_phi[phi_idx]/np.pi:.2f}pi")
        plt.savefig("%s/%s%s%05d_phi=" % (savedir, dname, aname, fnum) +f"{aa.possible_phi[phi_idx]/np.pi:.2f}pi"+".png")
    plt.close()

def tidal_profile_cmap(dname, fnum, grid_type):
    aname = "_tidal" #a for analysis
    tidal_profile.aname = aname
    data_location = file.data_loc + dname
    grid_type=grid_type
    savedir = file.savedir + dname + "/" + dname + aname + "/tidal_cmap"
    mkdir_if_not_exist(savedir)
    filename = "%s/disk.out1.%05d.athdf" % (data_location, fnum)
    aa = Athena_Analysis(filename=filename, grid_type=grid_type)
    
    aa.native_grid(get_r=True, get_z=True, get_phi=True)
    aa.get_primaries(get_rho=True, get_vel_phi=True)
    aa.get_potentials(get_accel=True, get_companion_grav=True)

    tidal_C = aa.native_to_cart(aa.get_C_vec(-1 * aa.rho * aa.gradient(aa.accel_pot + aa.companion_grav_pot, coordinates=grid_type)))
    lrl, lrl_cart = aa.get_lrl(components=True)
    eccent = aa.rho * lrl_cart
    tidal_growth = (vec.ortho_dot(eccent, tidal_C) / vec.get_magnitude(eccent)) / aa.integrate(aa.rho, "all")

    vert = 1
    horz = 2
    gs = gridspec.GridSpec(vert, horz)
    fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)
    ax_rho = fig.add_subplot(gs[0, 1])
    ax_tid = fig.add_subplot(gs[0, 0])
    
    aa.midplane_colorplot(tidal_growth, ax_tid, vbound=[-1e-4,1e-4], slicetype='z', log=False)
    ax_tid.set_title("Tidal Eccentricity Contribution")
    aa.midplane_colorplot(aa.rho, ax_rho, plot_COM=False, vbound=[1e-5, 1e2], slicetype='z')
    ax_rho.set_title(r"$\rho$")

    plt.tight_layout()
    plt.subplots_adjust(top=0.82)
    orbit = (aa.time / sim.binary_period)
   
    fig.suptitle(f"orbit: {orbit:.2f}")
    plt.savefig("%s/%s%s%05d.png" % (savedir, dname, aname, fnum))
    plt.close()

def reduce_dimensions(q, phi_idx, aa, mass_weight=False):
    if q.shape == aa.vector_array_size:
        return np.array([reduce_dimensions(q[0], phi_idx, aa), reduce_dimensions(q[1], phi_idx, aa), reduce_dimensions(q[2], phi_idx, aa)])
    if phi_idx == "az_ave":
        if not hasattr(reduce_dimensions, "normalization_weight"):
            if mass_weight:
                reduce_dimensions.normalization_weight = aa.integrate(aa.rho, "shell")
                print("this should print exactly once")
            else:
                reduce_dimensions.normalization_weight = aa.integrate(1, "shell")
                print("this should print exactly once")
        if mass_weight:
            return (aa.integrate(q*aa.rho, "shell") / reduce_dimensions.normalization_weight)
        else:
            return (aa.integrate(q, "shell") / reduce_dimensions.normalization_weight)
    else:
        if (not hasattr(reduce_dimensions, "normalization_weight")) and aa.gridtype == "Spherical":
            if mass_weight:
                reduce_dimensions.normalization_weight = aa.integrate(aa.rho, "theta")
                print("this should print exactly once")
            else:
                reduce_dimensions.normalization_weight = aa.integrate(1, "theta")
                print("this should print exactly once")
        if aa.gridtype == "Spherical":
            if mass_weight:
                return((aa.integrate(q*aa.rho, "theta") / reduce_dimensions.normalization_weight)[phi_idx])
            else:
                return((aa.integrate(q, "theta") / reduce_dimensions.normalization_weight)[phi_idx])
        if (not hasattr(reduce_dimensions, "normalization_weight")) and aa.gridtype == "Cylindrical":
            if mass_weight:
                reduce_dimensions.normalization_weight = aa.integrate(aa.rho, "z")
                print("this should print exactly once")
            else:
                reduce_dimensions.normalization_weight = aa.integrate(1, "z")
                print("this should print exactly once")
        if aa.gridtype == "Cylindrical":
            if mass_weight:
                return((aa.integrate(q*aa.rho, "z") / reduce_dimensions.normalization_weight)[phi_idx])
            else:
                return((aa.integrate(q, "z") / reduce_dimensions.normalization_weight)[phi_idx])
