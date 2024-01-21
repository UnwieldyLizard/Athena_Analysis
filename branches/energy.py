from .roots.athena_analysis import *
from .roots.misc_func import *
from datetime import datetime
import gc
import pickle

logging.basicConfig(filename=file.logs_loc+"/energy.log", encoding='utf-8', level=logging.INFO)

def kinetic_energy_plot(dnames, fnum_range, plot_every = 100, pickle_every = None):
    aname = "kinetic_energy" #a for analysis
    file_spacings = np.zeros(len(dnames), dtype=int)
    for d in range(len(dnames)):
        file_spacings[d] = find_file_spacing(dnames[d])
    savedir = file.savedir + "comparison" + "/" + aname
    pickldir = savedir + "/pickles"
    mkdir_if_not_exist(savedir)
    mkdir_if_not_exist(pickldir)
    file_spacing = np.gcd.reduce(file_spacings)
    now = datetime.now()

    if pickle_every is None:
        pickle_every = plot_every

    logging.info("setting up")
    found_load_point = False
    load_point = fnum_range[0]
    while found_load_point == False:
        if os.path.exists("%s/pickles/pickle_%05d.dat" % (savedir, load_point)):
            logging.info("Found data, loading from: %s" % load_point)
            found_load_point = True
        else:
            load_point -= file_spacing
        if load_point <= 0:
            load_point = 0
            logging.info("No load point found, restarting")
            break
    if load_point != 0:
        with open("%s/pickles/pickle_%05d.dat" % (savedir, load_point), "rb") as pickle_file:
            data = pickle.load(pickle_file)
        offset = len(data["orbits"])
    else:
        offset = 0

    orbits = np.zeros([(fnum_range[1] - fnum_range[0] + offset)])
    KEs = np.zeros([len(dnames), (fnum_range[1] - fnum_range[0] + offset), 3])
    KEBoundaries = np.zeros([len(dnames), (fnum_range[1] - fnum_range[0] + offset), 3])
    IntBoundaries = np.zeros([len(dnames), (fnum_range[1] - fnum_range[0] + offset), 3])
    final_points = np.zeros(len(dnames))

    if load_point != 0:
        if dnames != data["dnames"]:
            raise("this functionality doesn't exist yet, please use consistent dnames")
        orbits[:offset] = data["orbits"]
        KEs[:, :offset] = data["KEs"][:,:offset]
        KEBoundaries[:, :offset] = data["KEBoundaries"][:,:offset]
        IntBoundaries[:, :offset] = data["IntBoundaries"][:,:offset]
        final_points = data["final_points"]

    for i, f in enumerate(range(fnum_range[0], fnum_range[1], file_spacing)):
        logging.info(datetime.now()-now)
        now = datetime.now()
        logging.info("fnum = %d" % f)

        for d in range(len(dnames)):
            if final_points[d] == 0:
                logging.info(dnames[d])
                filename = "%s/disk.out1.%05d.athdf" % (file.data_loc + dnames[d], f)

                try: 
                    aa = Athena_Analysis(filename=filename, grid_type=file.grid_types[dnames[d]])
                except:
                    if final_points[d] == 0:
                        logging.info(f"{dnames[d]} ended at f = {f}")
                        print(f"{dnames[d]} ended at f = {f}")
                        final_points[d] = i
                    if np.all(final_points):
                        logging.info("All series terminated, pickling")
                        with open("%s/pickles/pickle_%05d.dat" % (savedir, f-1), "wb") as pickle_file:
                            pickle.dump({
                                "dnames": dnames,
                                "orbits": orbits[:i],
                                "KEs": KEs[:,:i],
                                "KEBoundaries": KEBoundaries[:,:i],
                                "IntBoundaries": IntBoundaries[:,:i],
                                "final_points": final_points,
                            }, pickle_file)
                        return
                    continue

                aa.get_primaries(get_rho=True, get_vel_r=True, get_vel_phi=True, get_vel_z=True)

                KEr = 0.5 * aa.rho * aa.vel_r * aa.vel_r
                KEphi = 0.5 * aa.rho * aa.vel_phi * aa.vel_phi
                KEz = 0.5 * aa.rho * aa.vel_z * aa.vel_z
                KEBoundaries[d, i, 0] = aa.get_boundary_flux(KEr)
                KEBoundaries[d, i, 1] = aa.get_boundary_flux(KEphi)
                KEBoundaries[d, i, 2] = aa.get_boundary_flux(KEz)
                IntBoundaries[d, i, 0] = np.sum(KEBoundaries[d,:i+1,0])
                IntBoundaries[d, i, 1] = np.sum(KEBoundaries[d,:i+1,1])
                IntBoundaries[d, i, 2] = np.sum(KEBoundaries[d,:i+1,2])
                KEs[d, i, 0] = aa.integrate(KEr, "All")
                KEs[d, i, 1] = aa.integrate(KEphi, "All")
                KEs[d, i, 2] = aa.integrate(KEz, "All")
                orbits[i] = aa.time / sim.binary_period

        if i % plot_every == 0:

            vert = 1
            horz = 3
            gs = gridspec.GridSpec(vert, horz)
            fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

            ax_KEr = fig.add_subplot(gs[0, 0])
            ax_KEr.set_title("Radial")
            ax_KEphi = fig.add_subplot(gs[0, 1])
            ax_KEphi.set_title("Azimuthal")
            ax_KEz = fig.add_subplot(gs[0, 2])
            ax_KEz.set_title("Vertical")

            for d in range(len(dnames)):
                if final_points[d] == 0:
                    ax_KEr.plot(orbits[:i+1], KEs[d, :i+1, 0] + IntBoundaries[d, :i+1, 0], f"C{d}-", label=dnames[d])
                    ax_KEphi.plot(orbits[:i+1], KEs[d, :i+1, 1] + IntBoundaries[d, :i+1, 1], f"C{d}-", label=dnames[d])
                    ax_KEz.plot(orbits[:i+1], KEs[d, :i+1, 2] + IntBoundaries[d, :i+1, 2], f"C{d}-", label=dnames[d])
                else:
                    ax_KEr.plot(orbits[:int(final_points[d])], KEs[d, :int(final_points[d]), 0] + IntBoundaries[d, :int(final_points[d]), 0], f"C{d}-", label=dnames[d])
                    ax_KEphi.plot(orbits[:int(final_points[d])], KEs[d, :int(final_points[d]), 1] + IntBoundaries[d, :int(final_points[d]), 1], f"C{d}-", label=dnames[d])
                    ax_KEz.plot(orbits[:int(final_points[d])], KEs[d, :int(final_points[d]), 2] + IntBoundaries[d, :int(final_points[d]), 2], f"C{d}-", label=dnames[d])
            
            plt.tight_layout()
            plt.legend()
            orbit = (aa.time / sim.binary_period)
            plt.subplots_adjust(top=(1-0.01*(16/vert)))
            fig.suptitle("Kinetic Energies Plus Cummulative KE Flux Through Boundaries")
            plt.savefig("%s/%s%05d.png" % (savedir, aname, f))
            plt.close()

        if i % pickle_every == 0:
            with open("%s/pickles/pickle_%05d.dat" % (savedir, f-1), "wb") as pickle_file:
                pickle.dump({
                    "dnames": dnames,
                    "orbits": orbits[:i],
                    "KEs": KEs[:,:i],
                    "KEBoundaries": KEBoundaries[:,:i],
                    "IntBoundaries": IntBoundaries[:,:i],
                    "final_points": final_points,
                }, pickle_file)
