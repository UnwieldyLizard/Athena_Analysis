from .roots.athena_analysis import *
from .roots.misc_func import *
import logging
import pickle

logging.basicConfig(filename=file.logs_loc+"/meccentricity_profile.log", encoding='utf-8', level=logging.INFO)


def eccentricity_profile(dname, fnum_limits, file_spacing=None, plot_every=100, pickle_every=None, pickle_flags=[], restart=False, alpha=None, stress_test = False):
    aname = "_eccent_growth_prec" #a for analysis
    if stress_test:
        aname += "_stress"
    data_location = file.data_loc + dname
    coordinates = file.grid_types[dname]
    MHD = file.MHD[dname]
    #data_location = file.mkitp + file.alpha3
    if file_spacing is None:
        file_spacing = find_file_spacing(dname)
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
        if stress_test:
            terms = ["tidal", "press", "boundary", "B", "vpress", "hpress", "vB", "hB"]
            sum_terms = ["tidal", "press", "boundary", "B"]
        else:
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
        
        if MAGNETIC_FIELDS_ENABLED:
            if stress_test:
                if aa.gridtype == "Sphereical":
                    logging.info("not implimented spherical stress calc in eccent prof")
                if aa.gridtype == "Cylindrical":
                    B_sq = ((aa.B_z * aa.B_z) + (aa.B_phi * aa.B_phi) + (aa.B_r * aa.B_r))
                    max_stress = np.array([[aa.B_z * aa.B_z - (B_sq/2), aa.B_z * aa.B_phi             , aa.B_z * aa.B_r],
                                            [aa.B_z * aa.B_phi        , aa.B_phi * aa.B_phi - (B_sq/2), aa.B_phi * aa.B_r],
                                            [aa.B_z * aa.B_r          , aa.B_r * aa.B_phi             , aa.B_r * aa.B_r - (B_sq/2)]])
                    force["B"] = aa.tensor_divergence(max_stress)
                    force["vB"] = np.array(force["B"])
                    force["vB"][1, :] *= 0
                    force["vB"][2, :] *= 0
                force["hB"] = force["B"] - force["vB"]
            else:    
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

        force["press"] = -1 * aa.gradient(aa.press, coordinates=coordinates)    
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

        #test_cart = aa.vec_vol_integrate(lrl_cart)

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

                if key in sum_terms:
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


def replot(dname, fnum, pname="", aspect_ratio=2, stress_test=False, recompute_sum=False):
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
    if stress_test:
        aname += "_stress"
    savedir = file.savedir + dname + "/" + dname + aname
    MHD = file.MHD[dname]
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
        if stress_test:
            terms = ["tidal", "press", "boundary", "B", "vpress", "hpress", "vB", "hB"]
            sum_terms = ["tidal", "press", "boundary", "B"]
        else:
            terms = ["tidal", "press", "boundary", "Bpress", "Btens", "vpress", "hpress", "vBpress", "hBpress"]
            sum_terms = ["tidal", "press", "boundary", "Bpress", "Btens"]
        MAGNETIC_FIELDS_ENABLED = True
    else:
        terms = ["tidal", "press", "boundary", "visc", "vpress", "hpress", "tidal_vr"]
        sum_terms = ["tidal", "press", "boundary", "visc"]
        MAGNETIC_FIELDS_ENABLED = False
    
    if recompute_sum:
        cutoff = len(data["integrated_eccent_sum_term_series"])
        integrated_eccent_sum_term_series = np.zeros(cutoff)
        integrated_eccent_phase_sum_term_series = np.zeros(cutoff)
        for key in sum_terms:
            integrated_eccent_sum_term_series += data["integrated_eccent_term_series"][key][:cutoff]
            integrated_eccent_phase_sum_term_series += data["integrated_eccent_term_phase_series"][key][:cutoff]

        initial_offset_diff = data["eccent_phase_series"][5] - integrated_eccent_phase_sum_term_series[5]
        integrated_eccent_phase_sum_term_series += initial_offset_diff
    else:
        integrated_eccent_sum_term_series = data["integrated_eccent_sum_term_series"]

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
    ax.plot(data["orbit_series"], integrated_eccent_sum_term_series, "C9--", label="total source terms")
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
    ax.plot(data["orbit_series"], (wrap_phase(integrated_eccent_phase_sum_term_series)) / np.pi, "C1--", label="total source terms")
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
    #ax.plot(data["orbit_series"], (wrap_phase(data["integrated_phase_flips"])) / np.pi, f"C3--", label="phase inversion")
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

def eccentricity_radial_profile(dname, fnum):
    aname = "_eccent_growth_prec" #a for analysis
    sname = "radial"
    data_location = file.data_loc + dname
    coordinates = file.grid_types[dname]
    MHD = file.MHD[dname]
    alpha = file.alpha[dname]
    #data_location = file.mkitp + file.alpha3
    savedir = file.savedir + dname + "/" + dname + aname
    mkdir_if_not_exist(savedir)
    mkdir_if_not_exist(savedir + "/pickles")

    if MHD == True:
        terms = ["tidal", "press", "Bpress", "Btens", "vpress", "hpress", "vBpress", "hBpress", "flux"]
        sum_terms = ["tidal", "press", "Bpress", "Btens", "flux"]
        MAGNETIC_FIELDS_ENABLED = True
    else:
        terms = ["tidal", "press", "visc", "vpress", "hpress", "tidal_vr", "flux"]
        sum_terms = ["tidal", "press", "visc", "flux"]
        MAGNETIC_FIELDS_ENABLED = False
    
    eccent_term = {}
    eccent_term_phase = {}

    filename = "%s/disk.out1.%05d.athdf" % (data_location, fnum)
    
    aa = Athena_Analysis(filename=filename, grid_type= coordinates)
    
    orbit = (aa.time / sim.binary_period)

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
    shell_mass = aa.integrate(aa.rho, variable='shell')
    rhoxlrl = np.array([aa.rho*lrl_cart[0], aa.rho*lrl_cart[1], np.zeros(aa.array_size)])
    shell_rhoxlrl = aa.vec_shell_integrate(rhoxlrl)
    mass_weighted_eccent = shell_rhoxlrl / shell_mass

    vel = aa.native_to_cart(np.array([aa.vel_z, aa.vel_phi, aa.vel_r]))
    eccent_flux = aa.rho * np.array([lrl_cart[0] * vel, lrl_cart[1] * vel, lrl_cart[2] * vel])  
    flux_term = aa.vec_shell_integrate(-1 * np.array([aa.divergence(eccent_flux[0]), aa.divergence(eccent_flux[1]), aa.divergence(eccent_flux[2])])) / shell_mass
    growth, prec = vec.growth_prec(flux_term, mass_weighted_eccent)
    eccent_term["flux"] = growth
    eccent_term_phase["flux"] = prec

    eccent = vec.get_magnitude(mass_weighted_eccent)
    eccent_phase = np.arctan2(mass_weighted_eccent[1], mass_weighted_eccent[0]) + sim.orbital_Omega * aa.time

    # calculate C
    aa.get_primaries(get_vel_r=True, get_vel_phi=True)
    C = {}
    for key in force:
        C[key] = aa.get_C_vec(force[key])
        C[key] = aa.native_to_cart(C[key])
        evol_C = aa.vec_shell_integrate(C[key]) / shell_mass
        growth, prec = vec.growth_prec(evol_C, mass_weighted_eccent)

        eccent_term[key] = growth
        eccent_term_phase[key] = prec

    flux = aa.integrate(aa.radial_transport(aa.rho * vec.get_magnitude(lrl_cart)), "Shell")

    aa.get_primaries(get_vel_r=True)
    aa.get_face_areas(get_rcc_face_areas=True)
    flux2 = aa.integrate(aa.rcc_face_area * aa.vel_r * aa.rho * vec.get_magnitude(lrl_cart), "Shell")

    # set up plot

    r_axis = aa.possible_r

    logging.info("\tplot")

    # growth plot
    vert = 3
    horz = 1
    gs = gridspec.GridSpec(vert, horz)
    fig = plt.figure(figsize=(2*horz*3, vert*3), dpi=300)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(r_axis, eccent, "C2-", label="eccent")
    ax.set_xlabel("r")
    ax.set_ylabel("eccent magnitude")
    #ax.set_ylim([1e-4, 1])
    #ax.set_yscale("log")
    #ax.set_title("mhd_3d eccent magnitude")
    ax.legend(loc="upper left")
    #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    ax = fig.add_subplot(gs[1, 0])
    for k, key in enumerate(terms):
        ax.plot(r_axis, eccent_term[key], f"C{k}-", label=key)
    ax.set_xlabel("r")
    ax.set_ylabel("eccent contribution")
    ax.set_ylim(top=25)
    #ax.set_title("mhd_3d time integrated sources")
    ax.legend(loc="upper left")
    #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    ax = fig.add_subplot(gs[2, 0])
    #ax.plot(r_axis, flux)
    ax.plot(r_axis, flux2)
    #ax.legend()
    ax.set_xlabel("r")
    ax.set_ylabel("Flux")
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    plt.subplots_adjust(top=(1-0.01*(16/vert)))
    title = f"{dname}: orbit {orbit:.2f}" 
    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig("%s/%s_growth_%05d_rad.png" % (savedir, dname, fnum))
    plt.close()

    # precession plot
    vert = 2
    horz = 1
    gs = gridspec.GridSpec(vert, horz)
    fig = plt.figure(figsize=(2*horz*3, vert*3), dpi=300)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(r_axis, wrap_phase(eccent_phase) / np.pi, "C3-", label="measured")
    ax.set_xlabel("r")
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
        ax.plot(r_axis, (wrap_phase(eccent_term_phase[key])) / np.pi, f"C{k}-", label=key)
    ax.set_xlabel("r")
    ax.set_ylabel("eccent phase contribution")
    ax.set_ylim([-1, 1])
    ax.set_title("mhd_3d time integrated sources")
    ax.legend(loc="upper left")
    #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    plt.subplots_adjust(top=(1-0.01*(16/vert)))
    title = f"{dname}: orbit {orbit:.2f}" 
    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig("%s/%s_prec_%05d_rad.png" % (savedir, dname, fnum))
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

def split_profile(dname, fnum_limits, radial_cutoff=5, plot_every=100, pickle_every=None, pickle_flags=[], restart=False, alpha=None):
    aname = "_eccent_growth_prec" #a for analysis
    sname = "_split"
    data_location = file.data_loc + dname
    coordinates = file.grid_types[dname]
    MHD = file.MHD[dname]
    #data_location = file.mkitp + file.alpha3
    savedir = file.savedir + dname + "/" + dname + aname
    mkdir_if_not_exist(savedir)
    mkdir_if_not_exist(savedir + "/"+sname+"_pickles")
    if pickle_every is None:
        pickle_every = plot_every
    file_spacing = find_file_spacing(dname)
    fnum_range = list(range(fnum_limits[0], fnum_limits[1], file_spacing))
    pickle_flags = np.array(pickle_flags)
    pickle_flags = np.append(pickle_flags, fnum_limits[1])
    if restart == True:
        dname = dname + "_%ssrt" % fnum_range[0]

    if MHD == True:
        terms = ["tidal", "press", "boundary", "cutoff flux", "B", "vpress", "hpress", "vB", "hB"]
        sum_terms = ["tidal", "press", "boundary", "cutoff flux", "B"]
        MAGNETIC_FIELDS_ENABLED = True
    else:
        terms = ["tidal", "press", "visc", "vpress", "hpress", "tidal_vr", "flux", "cutoff flux"]
        sum_terms = ["tidal", "press", "visc", "flux", "cutoff flux"]
        MAGNETIC_FIELDS_ENABLED = False
    
    orbit_series = np.zeros(len(fnum_range))
    time_series = np.zeros(len(fnum_range))
    in_eccent_series = np.zeros(len(fnum_range))
    in_eccent_phase_series = np.zeros(len(fnum_range))
    out_eccent_series = np.zeros(len(fnum_range))
    out_eccent_phase_series = np.zeros(len(fnum_range))

    in_eccent_term_series = {}
    in_integrated_eccent_term_series = {}
    in_eccent_term_phase_series = {}
    in_integrated_eccent_term_phase_series = {}
    out_eccent_term_series = {}
    out_integrated_eccent_term_series = {}
    out_eccent_term_phase_series = {}
    out_integrated_eccent_term_phase_series = {}

    for key in terms:
        in_eccent_term_series[key] = np.zeros(len(fnum_range))
        in_integrated_eccent_term_series[key] = np.zeros(len(fnum_range))
        in_eccent_term_phase_series[key] = np.zeros(len(fnum_range))
        in_integrated_eccent_term_phase_series[key] = np.zeros(len(fnum_range))
        out_eccent_term_series[key] = np.zeros(len(fnum_range))
        out_integrated_eccent_term_series[key] = np.zeros(len(fnum_range))
        out_eccent_term_phase_series[key] = np.zeros(len(fnum_range))
        out_integrated_eccent_term_phase_series[key] = np.zeros(len(fnum_range))
    in_integrated_phase_flips = np.zeros(len(fnum_range))
    out_integrated_phase_flips = np.zeros(len(fnum_range))

    in_integrated_eccent_sum_term_series = np.zeros(len(fnum_range))
    in_integrated_eccent_phase_sum_term_series = np.zeros(len(fnum_range))
    out_integrated_eccent_sum_term_series = np.zeros(len(fnum_range))
    out_integrated_eccent_phase_sum_term_series = np.zeros(len(fnum_range))

    now = datetime.now()
    for i, fnum in enumerate(fnum_range):

        if (i == 0) and (fnum != 0) and (restart == False):
            logging.info("setting up")
            found_load_point = False
            load_point = fnum - file_spacing
            while found_load_point == False:
                if os.path.exists("%s/%spickles/%s_pickle_%05d.dat" % (savedir, sname+"_", dname, load_point)):
                    logging.info("Found data, loading from: %s" % load_point)
                    found_load_point = True
                else:
                    load_point -= file_spacing
                if load_point <= 0:
                    raise("No load point, you need to restart")
            with open("%s/%spickles/%s_pickle_%05d.dat" % (savedir, sname+"_", dname, load_point), "rb") as pickle_file:
                data = pickle.load(pickle_file)
            i_0 = data["offset"]

            orbit_series = np.zeros(len(fnum_range)+i_0)
            time_series = np.zeros(len(fnum_range)+i_0)
            in_eccent_series = np.zeros(len(fnum_range)+i_0)
            in_eccent_phase_series = np.zeros(len(fnum_range)+i_0)
            out_eccent_series = np.zeros(len(fnum_range)+i_0)
            out_eccent_phase_series = np.zeros(len(fnum_range)+i_0)

            logging.info(str(i_0))
            logging.info(str(len(fnum_range)))

            for key in terms:
                in_eccent_term_series[key] = np.zeros(len(fnum_range)+i_0)
                in_integrated_eccent_term_series[key] = np.zeros(len(fnum_range)+i_0)
                in_eccent_term_phase_series[key] = np.zeros(len(fnum_range)+i_0)
                in_integrated_eccent_term_phase_series[key] = np.zeros(len(fnum_range)+i_0)
                out_eccent_term_series[key] = np.zeros(len(fnum_range)+i_0)
                out_integrated_eccent_term_series[key] = np.zeros(len(fnum_range)+i_0)
                out_eccent_term_phase_series[key] = np.zeros(len(fnum_range)+i_0)
                out_integrated_eccent_term_phase_series[key] = np.zeros(len(fnum_range)+i_0)
            in_integrated_phase_flips = np.zeros(len(fnum_range)+i_0)
            out_integrated_phase_flips = np.zeros(len(fnum_range)+i_0)

            in_integrated_eccent_sum_term_series = np.zeros(len(fnum_range)+i_0)
            in_integrated_eccent_phase_sum_term_series = np.zeros(len(fnum_range)+i_0)
            out_integrated_eccent_sum_term_series = np.zeros(len(fnum_range)+i_0)
            out_integrated_eccent_phase_sum_term_series = np.zeros(len(fnum_range)+i_0)

            time_series[:i_0+1] = data["time_series"]
            orbit_series[:i_0+1] = data["orbit_series"]
            in_eccent_series[:i_0+1] = data["in_eccent_series"]
            in_eccent_phase_series[:i_0+1] = data["in_eccent_phase_series"]
            out_eccent_series[:i_0+1] = data["out_eccent_series"]
            out_eccent_phase_series[:i_0+1] = data["out_eccent_phase_series"]
            for key in terms:
                in_eccent_term_series[key][:i_0+1] = data["in_eccent_term_series"][key][:i_0+1]
                in_eccent_term_phase_series[key][:i_0+1] = data["in_eccent_term_phase_series"][key][:i_0+1]
                in_integrated_eccent_term_series[key][:i_0+1] = data["in_integrated_eccent_term_series"][key][:i_0+1]
                in_integrated_eccent_term_phase_series[key][:i_0+1] = data["in_integrated_eccent_term_phase_series"][key][:i_0+1]
                out_eccent_term_series[key][:i_0+1] = data["out_eccent_term_series"][key][:i_0+1]
                out_eccent_term_phase_series[key][:i_0+1] = data["out_eccent_term_phase_series"][key][:i_0+1]
                out_integrated_eccent_term_series[key][:i_0+1] = data["out_integrated_eccent_term_series"][key][:i_0+1]
                out_integrated_eccent_term_phase_series[key][:i_0+1] = data["out_integrated_eccent_term_phase_series"][key][:i_0+1]
            in_integrated_eccent_sum_term_series[:i_0+1] = data["in_integrated_eccent_sum_term_series"]
            in_integrated_eccent_phase_sum_term_series[:i_0+1] = data["in_integrated_eccent_phase_sum_term_series"]
            in_integrated_phase_flips[:i_0+1] = data["in_integrated_phase_flips"][:i_0+1]
            out_integrated_eccent_sum_term_series[:i_0+1] = data["out_integrated_eccent_sum_term_series"]
            out_integrated_eccent_phase_sum_term_series[:i_0+1] = data["out_integrated_eccent_phase_sum_term_series"]
            out_integrated_phase_flips[:i_0+1] = data["out_integrated_phase_flips"][:i_0+1]


            in_old_mass_weighted_eccent = data["in_old_mass_weighted_eccent"]
            out_old_mass_weighted_eccent = data["out_old_mass_weighted_eccent"]
            
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

        #Mass Weighted Eccentricity
        aa.get_primaries(get_rho=True, get_press=True)
        aa.get_potentials(get_companion_grav=True, get_accel=True)
        aa.build_vectors()
        if MAGNETIC_FIELDS_ENABLED:
            aa.get_Bfields()
        force = {}

        force["tidal"] = -1 * aa.rho * aa.gradient(aa.accel_pot + aa.companion_grav_pot, coordinates=coordinates)
        
        if MAGNETIC_FIELDS_ENABLED:
            if aa.gridtype == "Sphereical":
                logging.info("not implimented spherical stress calc in eccent prof")
            if aa.gridtype == "Cylindrical":
                B_sq = ((aa.B_z * aa.B_z) + (aa.B_phi * aa.B_phi) + (aa.B_r * aa.B_r))
                max_stress = np.array([[aa.B_z * aa.B_z - (B_sq/2), aa.B_z * aa.B_phi             , aa.B_z * aa.B_r],
                                        [aa.B_z * aa.B_phi        , aa.B_phi * aa.B_phi - (B_sq/2), aa.B_phi * aa.B_r],
                                        [aa.B_z * aa.B_r          , aa.B_r * aa.B_phi             , aa.B_r * aa.B_r - (B_sq/2)]])
                force["B"] = aa.tensor_divergence(max_stress)
                force["vB"] = np.array(force["B"])
                force["vB"][1, :] *= 0
                force["vB"][2, :] *= 0
            force["hB"] = force["B"] - force["vB"]
        else:
            force["visc"] = aa.alpha_visc(alpha)
            logging.info("alpha visc in testing")

        force["press"] = -1 * aa.gradient(aa.press, coordinates=coordinates)    
        force["vpress"] = np.array(force["press"], dtype=np.float32)
        if aa.gridtype == "Spherical":
            force["vpress"][0, :] *= 0
            force["vpress"][2, :] *= 0
        elif aa.gridtype == "Cylindrical":
            force["vpress"][1, :] *= 0
            force["vpress"][2, :] *= 0
        force["hpress"] = force["press"] - force["vpress"]

        if i == 0:
            aa.native_grid(get_r=True)
            r_max = max(aa.possible_r)+1
            r_min = min(aa.possible_r)-1
            #overshoots endpoints by 1

        lrl_native, lrl_cart = aa.get_lrl(components = True)        
        in_mass = aa.integrate(aa.rho, variable='phi', second_variable='z', third_bounds=[r_min, radial_cutoff])
        out_mass = aa.integrate(aa.rho, variable='phi', second_variable='z', third_bounds=[radial_cutoff, r_max])
        rhoxlrl = np.array([aa.rho*lrl_cart[0], aa.rho*lrl_cart[1], np.zeros(aa.array_size)])
        in_rhoxlrl = aa.vec_vol_integrate(rhoxlrl, radial_bounds=[r_min, radial_cutoff])
        out_rhoxlrl = aa.vec_vol_integrate(rhoxlrl, radial_bounds=[radial_cutoff, r_max])
        in_mass_weighted_eccent = in_rhoxlrl / in_mass
        out_mass_weighted_eccent = out_rhoxlrl / out_mass

        in_eccent_series[j] = vec.get_magnitude(in_mass_weighted_eccent)
        in_eccent_phase_series[j] = np.arctan2(in_mass_weighted_eccent[1], in_mass_weighted_eccent[0]) + sim.orbital_Omega * aa.time
        out_eccent_series[j] = vec.get_magnitude(out_mass_weighted_eccent)
        out_eccent_phase_series[j] = np.arctan2(out_mass_weighted_eccent[1], out_mass_weighted_eccent[0]) + sim.orbital_Omega * aa.time
        
        # calculate C
        aa.get_primaries(get_vel_r=True, get_vel_phi=True)
        C = {}
        for key in force:
            C[key] = aa.get_C_vec(force[key])
            C[key] = aa.native_to_cart(C[key])
            in_evol_C = aa.vec_vol_integrate(C[key], radial_bounds=[r_min, radial_cutoff]) / in_mass
            out_evol_C = aa.vec_vol_integrate(C[key], radial_bounds=[radial_cutoff, r_max]) / out_mass
            
            growth, prec = vec.growth_prec(in_evol_C, in_mass_weighted_eccent)
            in_eccent_term_series[key][j] = growth
            in_eccent_term_phase_series[key][j] = prec

            growth, prec = vec.growth_prec(out_evol_C, out_mass_weighted_eccent)
            out_eccent_term_series[key][j] = growth
            out_eccent_term_phase_series[key][j] = prec

        #Boundary
        tot_flux_0, in_flux_0, out_flux_0 = aa.get_boundary_flux(rhoxlrl[0], intermediates=True)
        tot_flux_1, in_flux_1, out_flux_1 = aa.get_boundary_flux(rhoxlrl[1], intermediates=True)
        rhoxlrl_flux_outer_boundary = np.array([out_flux_0, out_flux_1, 0])
        rhoxlrl_flux_inner_boundary = np.array([in_flux_0, in_flux_1, 0])
        advect_evol_outer_boundary = -1 * rhoxlrl_flux_outer_boundary / out_mass
        advect_evol_inner_boundary = -1 * rhoxlrl_flux_inner_boundary / in_mass
        adv_growth_outer_boundary, adv_prec_outer_boundary = vec.growth_prec(advect_evol_outer_boundary, out_mass_weighted_eccent)
        adv_growth_inner_boundary, adv_prec_inner_boundary = vec.growth_prec(advect_evol_inner_boundary, in_mass_weighted_eccent)

        mass_flux_boundary, inner_mass_flux, outer_mass_flux = aa.get_boundary_flux(aa.rho, intermediates=True)

        massd_evol_outer = outer_mass_flux * np.array(out_mass_weighted_eccent) / out_mass
        massd_evol_inner = inner_mass_flux * np.array(in_mass_weighted_eccent) / in_mass
        
        growth, prec = vec.growth_prec(massd_evol_inner, in_mass_weighted_eccent)
        in_eccent_term_series["boundary"][j] = adv_growth_inner_boundary + growth
        in_eccent_term_phase_series["boundary"][j] = adv_prec_inner_boundary + prec

        growth, prec = vec.growth_prec(massd_evol_outer, out_mass_weighted_eccent)
        out_eccent_term_series["boundary"][j] = adv_growth_outer_boundary + growth
        out_eccent_term_phase_series["boundary"][j] = adv_prec_outer_boundary + prec

        #Cutoff Boundary
        rhoxlrl_flux_cutoff = np.array([aa.get_slice_flux(rhoxlrl[0], radial_cutoff), aa.get_slice_flux(rhoxlrl[1], radial_cutoff), 0])
        out_advect_evol_cutoff = -1 * rhoxlrl_flux_cutoff / out_mass
        in_advect_evol_cutoff = +1 * rhoxlrl_flux_cutoff / in_mass
        out_adv_growth_cutoff, out_adv_prec_cutoff = vec.growth_prec(out_advect_evol_cutoff, out_mass_weighted_eccent)
        in_adv_growth_cutoff, in_adv_prec_cutoff = vec.growth_prec(in_advect_evol_cutoff, in_mass_weighted_eccent)
        
        mass_flux = aa.get_slice_flux(aa.rho, radial_cutoff)

        out_massd_evol_cutoff = -1 * mass_flux * np.array(out_mass_weighted_eccent) / out_mass
        in_massd_evol_cutoff = mass_flux * np.array(in_mass_weighted_eccent) / in_mass
        
        growth, prec = vec.growth_prec(in_massd_evol_cutoff, in_mass_weighted_eccent)
        in_eccent_term_series["cutoff flux"][j] = in_adv_growth_cutoff + growth
        in_eccent_term_phase_series["cutoff flux"][j] = in_adv_prec_cutoff + prec

        growth, prec = vec.growth_prec(out_massd_evol_cutoff, out_mass_weighted_eccent)
        out_eccent_term_series["cutoff flux"][j] = out_adv_growth_cutoff + growth
        out_eccent_term_phase_series["cutoff flux"][j] = out_adv_prec_cutoff + prec

        #integrate values forward
        in_total_eccent_increment = 0
        in_total_phase_increment = 0
        out_total_eccent_increment = 0
        out_total_phase_increment = 0
        for key in terms:
            if j != 0:
                eccent_increment = (time_series[j]-time_series[j-1])*(in_eccent_term_series[key][j] + in_eccent_term_series[key][j-1])/2
                in_integrated_eccent_term_series[key][j] = in_integrated_eccent_term_series[key][j-1] + eccent_increment
                phase_increment = (time_series[j]-time_series[j-1])*(in_eccent_term_phase_series[key][j] + in_eccent_term_phase_series[key][j-1])/2
                in_integrated_eccent_term_phase_series[key][j] = in_integrated_eccent_term_phase_series[key][j-1] + phase_increment

                if key in sum_terms:
                    in_total_eccent_increment += eccent_increment
                    in_total_phase_increment += phase_increment

                eccent_increment = (time_series[j]-time_series[j-1])*(out_eccent_term_series[key][j] + out_eccent_term_series[key][j-1])/2
                out_integrated_eccent_term_series[key][j] = out_integrated_eccent_term_series[key][j-1] + eccent_increment
                phase_increment = (time_series[j]-time_series[j-1])*(out_eccent_term_phase_series[key][j] + out_eccent_term_phase_series[key][j-1])/2
                out_integrated_eccent_term_phase_series[key][j] = out_integrated_eccent_term_phase_series[key][j-1] + phase_increment

                if key in sum_terms:
                    out_total_eccent_increment += eccent_increment
                    out_total_phase_increment += phase_increment

        if j == 0 or j == 1:
            in_integrated_eccent_sum_term_series[i] = in_eccent_series[i]
            in_integrated_eccent_phase_sum_term_series[i] = in_eccent_phase_series[i]
            out_integrated_eccent_sum_term_series[i] = out_eccent_series[i]
            out_integrated_eccent_phase_sum_term_series[i] = out_eccent_phase_series[i]
        elif j != 0 and j != 1:
            in_integrated_eccent_sum_term_series[j] = in_integrated_eccent_sum_term_series[j-1] + in_total_eccent_increment
            out_integrated_eccent_sum_term_series[j] = out_integrated_eccent_sum_term_series[j-1] + out_total_eccent_increment
            
            #checking for flip and adjusting
            if in_integrated_eccent_sum_term_series[j] < 0:
                in_flip_phase = np.pi
                in_integrated_eccent_sum_term_series[j] = np.abs(in_integrated_eccent_sum_term_series[j])
                logging.info("auto-flipped phase")
            elif (vec.ortho_dot(in_old_mass_weighted_eccent, in_mass_weighted_eccent)) < 0:
                in_flip_phase = np.pi
                logging.info("manually flipped phase")
            else:
                in_flip_phase = 0
            in_integrated_phase_flips[j] = in_integrated_phase_flips[j-1] + in_flip_phase
            in_integrated_eccent_phase_sum_term_series[j] = in_integrated_eccent_phase_sum_term_series[j-1] + in_total_phase_increment + in_flip_phase

            if out_integrated_eccent_sum_term_series[j] < 0:
                out_flip_phase = np.pi
                out_integrated_eccent_sum_term_series[j] = np.abs(out_integrated_eccent_sum_term_series[j])
                logging.info("auto-flipped phase")
            elif (vec.ortho_dot(out_old_mass_weighted_eccent, out_mass_weighted_eccent)) < 0:
                out_flip_phase = np.pi
                logging.info("manually flipped phase")
            else:
                out_flip_phase = 0
            out_integrated_phase_flips[j] = out_integrated_phase_flips[j-1] + out_flip_phase
            out_integrated_eccent_phase_sum_term_series[j] = out_integrated_eccent_phase_sum_term_series[j-1] + out_total_phase_increment + out_flip_phase
        
        #memorising old eccentricity vector for comparison in next time step
        in_old_mass_weighted_eccent = in_mass_weighted_eccent
        out_old_mass_weighted_eccent = out_mass_weighted_eccent

        # set up plot
        if j % plot_every == 0 or i == 0:
            logging.info("\tplot")

            # growth plot
            vert = 2
            horz = 2
            gs = gridspec.GridSpec(vert, horz)
            fig = plt.figure(figsize=(2*horz*3, vert*3), dpi=300)

            ax = fig.add_subplot(gs[0, 0])
            ax.set_title(f"Inner Disk (inside r = {10})")
            ax.plot(orbit_series[:j+1], in_eccent_series[:j+1], "C2-", label="measured")
            ax.plot(orbit_series[:j+1], in_integrated_eccent_sum_term_series[:j+1], "C9--", label="sum terms")
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
                ax.plot(orbit_series[:j+1], in_integrated_eccent_term_series[key][:j+1], f"C{k}-", label=key)
            ax.set_xlabel("binary orbit")
            ax.set_ylabel("eccent contribution")
            #ax.set_title("mhd_3d time integrated sources")
            ax.legend(loc="upper left")
            #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

            ax = fig.add_subplot(gs[0, 1])
            ax.set_title(f"Outer Disk (outside r = {10})")
            ax.plot(orbit_series[:j+1], out_eccent_series[:j+1], "C2-", label="measured")
            ax.plot(orbit_series[:j+1], out_integrated_eccent_sum_term_series[:j+1], "C9--", label="sum terms")
            ax.set_xlabel("binary orbit")
            ax.set_ylabel("eccent magnitude")
            ax.set_ylim([1e-4, 1])
            ax.set_yscale("log")
            #ax.set_title("mhd_3d eccent magnitude")
            ax.legend(loc="upper left")
            #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

            ax = fig.add_subplot(gs[1, 1])
            for k, key in enumerate(terms):
                ax.plot(orbit_series[:j+1], out_integrated_eccent_term_series[key][:j+1], f"C{k}-", label=key)
            ax.set_xlabel("binary orbit")
            ax.set_ylabel("eccent contribution")
            #ax.set_title("mhd_3d time integrated sources")
            ax.legend(loc="upper left")
            #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

            plt.tight_layout()
            plt.savefig("%s/%s_growth_%05d%s.png" % (savedir, dname, fnum, sname))
            plt.close()

            # precession plot
            vert = 2
            horz = 2
            gs = gridspec.GridSpec(vert, horz)
            fig = plt.figure(figsize=(2*horz*3, vert*3), dpi=300)

            ax = fig.add_subplot(gs[0, 0])
            ax.set_title(f"Inner Disk (inside r = {10})")
            ax.plot(orbit_series[:j+1], wrap_phase(in_eccent_phase_series[:j+1]) / np.pi, "C3-", label="measured")
            ax.plot(orbit_series[:j+1], (wrap_phase(in_integrated_eccent_phase_sum_term_series[:j+1])) / np.pi, "C1--", label="total source terms")
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
                ax.plot(orbit_series[:j+1], (wrap_phase(in_integrated_eccent_term_phase_series[key][:j+1])) / np.pi, f"C{k}-", label=key)
            ax.plot(orbit_series[:j+1], (wrap_phase(in_integrated_phase_flips[:j+1])) / np.pi, f"C3--", label="phase inversion")
            ax.set_xlabel("binary orbit")
            ax.set_ylabel("eccent phase contribution")
            ax.set_ylim([-1, 1])
            ax.set_title("mhd_3d time integrated sources")
            ax.legend(loc="upper left")
            #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
            #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

            ax = fig.add_subplot(gs[0, 1])
            ax.set_title(f"Outer Disk (outside r = {10})")
            ax.plot(orbit_series[:j+1], wrap_phase(out_eccent_phase_series[:j+1]) / np.pi, "C3-", label="measured")
            ax.plot(orbit_series[:j+1], (wrap_phase(out_integrated_eccent_phase_sum_term_series[:j+1])) / np.pi, "C1--", label="total source terms")
            ax.set_xlabel("binary orbit")
            ax.set_ylabel("eccent phase")
            ax.set_ylim([-1, 1])
            #ax.set_yscale("log")
            ax.set_title("mhd_3d eccent phase")
            ax.legend(loc="upper left")
            #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
            #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

            ax = fig.add_subplot(gs[1, 1])
            for k, key in enumerate(terms):
                ax.plot(orbit_series[:j+1], (wrap_phase(out_integrated_eccent_term_phase_series[key][:j+1])) / np.pi, f"C{k}-", label=key)
            ax.plot(orbit_series[:j+1], (wrap_phase(out_integrated_phase_flips[:j+1])) / np.pi, f"C3--", label="phase inversion")
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
            with open("%s/%spickles/%s_pickle_%05d.dat" % (savedir, sname+"_", dname, fnum), "wb") as pickle_file:
                pickle.dump({
                    "offset": j,
                    "time_series": time_series[:j+1],
                    "orbit_series": orbit_series[:j+1],
                    "in_eccent_series": in_eccent_series[:j+1],
                    "in_eccent_phase_series": in_eccent_phase_series[:j+1],
                    "out_eccent_series": out_eccent_series[:j+1],
                    "out_eccent_phase_series": out_eccent_phase_series[:j+1],
                    "in_eccent_term_series": in_eccent_term_series,
                    "in_eccent_term_phase_series": in_eccent_term_phase_series,
                    "out_eccent_term_series": out_eccent_term_series,
                    "out_eccent_term_phase_series": out_eccent_term_phase_series,
                    "in_integrated_eccent_term_series": in_integrated_eccent_term_series,
                    "in_integrated_eccent_term_phase_series": in_integrated_eccent_term_phase_series,
                    "out_integrated_eccent_term_series": out_integrated_eccent_term_series,
                    "out_integrated_eccent_term_phase_series": out_integrated_eccent_term_phase_series,
                    "in_integrated_eccent_sum_term_series": in_integrated_eccent_sum_term_series[:j+1],
                    "in_integrated_eccent_phase_sum_term_series": in_integrated_eccent_phase_sum_term_series[:j+1],
                    "out_integrated_eccent_sum_term_series": out_integrated_eccent_sum_term_series[:j+1],
                    "out_integrated_eccent_phase_sum_term_series": out_integrated_eccent_phase_sum_term_series[:j+1],
                    "in_integrated_phase_flips": in_integrated_phase_flips[:j+1],
                    "in_old_mass_weighted_eccent": in_old_mass_weighted_eccent,
                    "out_integrated_phase_flips": out_integrated_phase_flips[:j+1],
                    "out_old_mass_weighted_eccent": out_old_mass_weighted_eccent,
                }, pickle_file)

def split_profile_replot(dname, fnum, pname="", aspect_ratio=2):
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
    sname = "_split"
    savedir = file.savedir + dname + "/" + dname + aname
    logging.info("looking for file")
    found_load_point = False
    load_point = fnum
    while found_load_point == False:
        if os.path.exists("%s/%spickles/%s_pickle_%05d.dat" % (savedir, sname+"_", dname, fnum)):
            logging.info("Found data, plotting from: %s" % load_point)
            found_load_point = True
        else:
            load_point -= 1
        if load_point <= 0:
            raise("No data found to plot")
    with open("%s/%spickles/%s_pickle_%05d.dat" % (savedir, sname+"_", dname, fnum), "rb") as pickle_file:
        data = pickle.load(pickle_file)

    # growth plot
    vert = 2
    horz = 1

    if aspect_ratio >= 1:
        vert_scale = 1
        horz_scale = aspect_ratio
    else:
        vert_scale = 1 / aspect_ratio
        horz_scale = 1

    gs = gridspec.GridSpec(vert, horz)
    fig = plt.figure(figsize=(horz_scale * horz*3, vert_scale * vert*3), dpi=300)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(data["orbit_series"], data["in_eccent_series"], "C2-", label="inner disk")
    ax.plot(data["orbit_series"], data["out_eccent_series"], "C9-", label="outer disk")
    ax.set_xlabel("binary orbit")
    ax.set_ylabel("eccent magnitude")
    ax.set_ylim([1e-4, 1])
    ax.set_yscale("log")
    #ax.set_title("mhd_3d eccent magnitude")
    ax.legend(loc="upper left")
    #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    ax = fig.add_subplot(gs[1, 0])
    ax.plot(data["orbit_series"], wrap_phase(data["in_eccent_phase_series"]) / np.pi, "C3-", label="inner disk")
    ax.plot(data["orbit_series"], wrap_phase(data["out_eccent_phase_series"]) / np.pi, "C1-", label="outer disk")
    ax.set_xlabel("binary orbit")
    ax.set_ylabel("eccent phase")
    ax.set_ylim([-1, 1])
    #ax.set_yscale("log")
    ax.set_title("mhd_3d eccent phase")
    ax.legend(loc="upper left")
    #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    plt.tight_layout()
    plt.savefig("%s/%s%s_growth_%05d%s.png" % (savedir, pname, dname, fnum, sname))
    plt.close()

def split_profile_rates(dname, fnum, pname="", aspect_ratio=2, num_orbits_ave=1):
    """
    replots data + precession rates

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
    sname = "_split"
    savedir = file.savedir + dname + "/" + dname + aname
    logging.info("looking for file")
    found_load_point = False
    load_point = fnum
    while found_load_point == False:
        if os.path.exists("%s/%spickles/%s_pickle_%05d.dat" % (savedir, sname+"_", dname, fnum)):
            logging.info("Found data, plotting from: %s" % load_point)
            found_load_point = True
        else:
            load_point -= 1
        if load_point <= 0:
            raise("No data found to plot")
    with open("%s/%spickles/%s_pickle_%05d.dat" % (savedir, sname+"_", dname, fnum), "rb") as pickle_file:
        data = pickle.load(pickle_file)

    num_files_ave = int(num_orbits_ave * sim.filenums_per_orbit)

    # growth plot
    vert = 2
    horz = 3

    if aspect_ratio >= 1:
        vert_scale = 1
        horz_scale = aspect_ratio
    else:
        vert_scale = 1 / aspect_ratio
        horz_scale = 1

    gs = gridspec.GridSpec(vert, horz)
    fig = plt.figure(figsize=(horz_scale * horz*3, vert_scale * vert*3), dpi=300)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(data["orbit_series"], data["out_eccent_series"], "C9-", label="outer disk")
    ax.plot(data["orbit_series"], data["in_eccent_series"], "C2-", label="inner disk")
    ax.set_xlabel("binary orbit")
    ax.set_ylabel("eccent magnitude")
    ax.set_ylim([1e-4, 1])
    ax.set_yscale("log")
    ax.set_title("eccent magnitude")
    ax.legend(loc="upper left")
    #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    ax = fig.add_subplot(gs[0, 1])
    dorbit = data["orbit_series"][1] - data["orbit_series"][0]
    deccent_in = data["in_eccent_series"][1:] - data["in_eccent_series"][:-1]
    deccent_out = data["out_eccent_series"][1:] - data["out_eccent_series"][:-1]
    growth_rate_in = deccent_in/dorbit
    growth_rate_out = deccent_out/dorbit
    rate_orbit_axis = data["orbit_series"][1:] - dorbit
    ax.plot(rate_orbit_axis, growth_rate_out, "C9-", label="outer disk")
    ax.plot(rate_orbit_axis, growth_rate_in, "C2-", label="inner disk")
    ax.set_xlabel("binary orbit")
    ax.set_ylabel("eccent growth rate")
    ax.set_ylim([-0.5, 0.5])
    #ax.set_yscale("log")
    ax.set_title("eccent magnitude growth rate")
    ax.legend(loc="upper left")
    #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g /orbit'))

    ax = fig.add_subplot(gs[0, 2])
    growth_rate_in_ave = (growth_rate_in[num_files_ave:])
    growth_rate_out_ave = (growth_rate_out[num_files_ave:])
    orbit_ave_axis = (rate_orbit_axis[num_files_ave:])
    for k in range(num_files_ave-1):
        growth_rate_in_ave = growth_rate_in_ave + (growth_rate_in[num_files_ave-(k+1):-(k+1)])
        growth_rate_out_ave = growth_rate_out_ave + (growth_rate_out[num_files_ave-(k+1):-(k+1)])
        orbit_ave_axis = orbit_ave_axis + (rate_orbit_axis[num_files_ave-(k+1):-(k+1)])
    growth_rate_in_ave = growth_rate_in_ave / num_files_ave
    growth_rate_out_ave = growth_rate_out_ave / num_files_ave
    orbit_ave_axis = orbit_ave_axis / num_files_ave
    ax.plot(orbit_ave_axis, growth_rate_out_ave, "C9-", label="outer disk")
    ax.plot(orbit_ave_axis, growth_rate_in_ave, "C2-", label="inner disk")
    ax.set_xlabel("binary orbit")
    ax.set_ylabel("eccent growth rate")
    ax.set_ylim([-0.05, 0.05])
    #ax.set_yscale("log")
    ax.set_title(f"eccent magnitude growth rate averaged over {num_orbits_ave} orbits")
    ax.legend(loc="upper left")
    #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g /orbit'))

    ax = fig.add_subplot(gs[1, 0])
    ax.plot(data["orbit_series"], wrap_phase(data["out_eccent_phase_series"]) / np.pi, "C1-", label="outer disk")
    ax.plot(data["orbit_series"], wrap_phase(data["in_eccent_phase_series"]) / np.pi, "C3-", label="inner disk")
    ax.set_xlabel("binary orbit")
    ax.set_ylabel("eccent phase")
    ax.set_ylim([-1, 1])
    #ax.set_yscale("log")
    ax.set_title("eccent phase")
    ax.legend(loc="upper left")
    #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$'))
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    ax = fig.add_subplot(gs[1, 1])
    deccent_ph_in = data["in_eccent_phase_series"][1:] - data["in_eccent_phase_series"][:-1]
    deccent_ph_out = data["out_eccent_phase_series"][1:] - data["out_eccent_phase_series"][:-1]
    flips_bool_in = (deccent_ph_in >= 1.9*np.pi)
    flips_bool_out = (deccent_ph_out >= 1.9*np.pi)
    neg_flips_bool_in = (deccent_ph_in <= -1.9*np.pi)
    neg_flips_bool_out = (deccent_ph_out <= -1.9*np.pi)
    deccent_ph_in = deccent_ph_in - 2*np.pi*flips_bool_in + 2*np.pi*neg_flips_bool_in
    deccent_ph_out = deccent_ph_out - 2*np.pi*flips_bool_out + 2*np.pi*neg_flips_bool_out
    prec_rate_in = deccent_ph_in#/dorbit
    prec_rate_out = deccent_ph_out#/dorbit
    ax.plot(rate_orbit_axis, prec_rate_out / np.pi, "C1-", label="outer disk")
    ax.plot(rate_orbit_axis, prec_rate_in / np.pi, "C3-", label="inner disk")
    ax.set_xlabel("binary orbit")
    ax.set_ylabel("eccent phase precession rate")
    ax.set_ylim([-0.1, 0.1])
    #ax.set_yscale("log")
    ax.set_title("eccent phase precession rate")
    ax.legend(loc="upper left")
    #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$/orbit'))
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    ax = fig.add_subplot(gs[1, 2])
    prec_rate_in_ave = (prec_rate_in[num_files_ave:])
    prec_rate_out_ave = (prec_rate_out[num_files_ave:])
    for k in range(num_files_ave-1):
        prec_rate_in_ave = prec_rate_in_ave + (prec_rate_in[num_files_ave-(k+1):-(k+1)])
        prec_rate_out_ave = prec_rate_out_ave + (prec_rate_out[num_files_ave-(k+1):-(k+1)])
    prec_rate_in_ave = prec_rate_in_ave / num_files_ave
    prec_rate_out_ave = prec_rate_out_ave / num_files_ave
    ax.plot(orbit_ave_axis, prec_rate_out_ave / np.pi, "C1-", label="outer disk")
    ax.plot(orbit_ave_axis, prec_rate_in_ave / np.pi, "C3-", label="inner disk")
    ax.set_xlabel("binary orbit")
    ax.set_ylabel("eccent phase precession rate")
    ax.set_ylim([-0.01, 0.01])
    #ax.set_yscale("log")
    ax.set_title(f"eccent phase precession rate averaged over {num_orbits_ave} orbits")
    ax.legend(loc="upper left")
    #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g $\pi$/orbit'))
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2g'))

    plt.subplots_adjust(top=(1-0.01*(16/vert)))
    plt.suptitle(dname)
    plt.tight_layout()
    plt.savefig("%s/%s%02d%s_growth_%05d%s.png" % (savedir, pname, num_orbits_ave, dname, fnum, sname))
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
