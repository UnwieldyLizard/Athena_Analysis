from roots.athena_analysis import *
from roots.misc_func import *
import logging
import scipy.fft

#logging.basicConfig(filename='/home/morgan/mresearch/testing.log', encoding='utf-8', level=logging.INFO)

def main(dname, fnum_range, file_spacing=1):
    #aname = "_alpha"
    aname = "_deviation"
    #data_location = "/home/morgan/mnt/kitp/data/cvdisk/superhump_3d_alpha03"
    #data_location = "/home/morgan/mnt/kitp/data2/cvdisk/CVThin2/Data"
    data_location = file.data_loc + dname
    #grid_type="Spherical"
    grid_type="Cylindrical"
    #data_location = "/home/morgan/mnt/kitp/data/cvdisk/CVThick3"
    savedir = file.savedir + dname + "/" + dname + aname
    mkdir_if_not_exist(savedir)

    for i, fnum in enumerate(range(fnum_range[0], fnum_range[1], file_spacing)):
        print("fnum = %d" % fnum)
        filename = "%s/disk.out1.%05d.athdf" % (data_location, fnum)

        aa = Athena_Analysis(filename=filename, grid_type=grid_type)

        #calculate stuff

        aa.get_primaries(get_rho=True)
        aa.native_grid(get_r=True)

        deviations = np.copy(aa.rho)
        means = (aa.integrate(aa.rho, "shell") / aa.integrate(1, "shell"))
        
        for n in range(aa.NumMeshBlocks):
            for r in range(aa.r_len):
                r_loc = int(np.searchsorted(aa.possible_r, aa.r_primitive[n, r]))
                deviations[n,:,:,r] -= means[r_loc]

        #plot stuff 
        vert = 1
        horz = 2
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        
        rotation = -1*sim.orbital_Omega * aa.time

        aa.midplane_colorplot(aa.rho, ax1, rotation=rotation)
        aa.midplane_colorplot(deviations, ax2, vbound=[-1,1], rotation=rotation)
        #ax1.plot(aa.possible_r, aa.integrate(aa.rho, "shell")/aa.integrate(1, "shell"))
       
        ax1.set_title(r"Density")
        ax2.set_title(r"Density Deviation From Shell Mean")
        #ax1.set_ylabel(r"$\rho$")
        #ax1.set_xlabel(r"r")
       

        plt.tight_layout()
        plt.subplots_adjust(top=(1-0.01*(16/vert)))
        orbit = (aa.time / sim.binary_period)
        fig.suptitle(f"Orbit: {orbit:.2f}")
        plt.savefig("%s/%s%s_test_plot%05d.png" % (savedir, dname, aname, fnum))
        plt.close()


main("Cyl_15_2", [4000, 4003])

def ff2_test():
    savedir = file.savedir + "test"
    field = np.zeros((100, 100))
    x = np.arange(0, 2, 2/100)
    for n in range(100):
        field[n] = np.exp((x*1j * np.pi))

    vert = 1
    horz = 3
    gs = gridspec.GridSpec(vert, horz)
    fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    im1 = ax1.imshow(field)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad=0.05)
    ax1.grid(False)
    cbar = plt.colorbar(im1, cax=cax1)
    cbar.ax.tick_params()

    fft2 = scipy.fft.fft2(field)

    im2 = ax2.imshow(np.abs(scipy.fft.fftshift(fft2)))
    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes("right", size="5%", pad=0.05)
    ax2.grid(False)
    cbar = plt.colorbar(im2, cax=cax2)
    cbar.ax.tick_params()

    im3 = ax3.imshow(np.real(scipy.fft.ifft2(fft2)))
    divider = make_axes_locatable(ax3)
    cax3 = divider.append_axes("right", size="5%", pad=0.05)
    ax3.grid(False)
    cbar = plt.colorbar(im3, cax=cax3)
    cbar.ax.tick_params()

    plt.tight_layout()
    plt.savefig("%s/ff2_test.png" % (savedir))

#ff2_test()