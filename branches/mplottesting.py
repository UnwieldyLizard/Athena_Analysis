from .roots.athena_analysis import *
import logging
import scipy.fft

logging.basicConfig(filename='/home/morgan/mresearch/testing.log', encoding='utf-8', level=logging.INFO)

def main():
    aname = "Press"
    dname = "Cyl_1"
    #data_location = "/home/morgan/mnt/kitp/data/cvdisk/superhump_3d_alpha03"
    #data_location = "/home/morgan/mnt/kitp/data2/cvdisk/CVThin2/Data"
    data_location = file.data_loc + dname
    #grid_type="Spherical"
    grid_type="Cylindrical"
    #data_location = "/home/morgan/mnt/kitp/data/cvdisk/CVThick3"
    savedir = file.savedir + "test"
    mkdir_if_not_exist(savedir)

    for i, fnum in enumerate(range(5000, 5004, 2)):
        print("fnum = %d" % fnum)
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

        if i != 0:
            aa.r = r_grid
            aa.phi = phi_grid
            aa.phi_len = phi_len
            if grid_type == "Spherical":
                aa.theta_len = theta_len
            if grid_type == "Cylindrical":
                aa.z_len = z_len
            aa.r_len = r_len
        
        #calculate stuff
        #aa.get_primaries(get_press=True, get_vel_r=True, get_vel_phi=True)
        #aa.get_potentials(get_companion_grav=True, get_accel=True)
        #tidal = -1 * aa.gradient(aa.companion_grav_pot + aa.accel_pot, coordinates='Spherical')
        aa.get_primaries(get_rho=True)
        #aa.native_grid(get_r=True, get_z=True, get_phi=True)
        #drhodr = aa.differentiate(aa.rho, 'r')
        #drhodtheta = aa.differentiate(aa.rho, 'theta')
        #drhodphi = aa.differentiate(aa.rho, 'phi')
        #drhodz = aa.differentiate(aa.rho, 'z')

        #plot stuff 
        vert = 2
        horz = 2
        gs = gridspec.GridSpec(vert, horz)
        fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

        ax1 = fig.add_subplot(gs[0, 1])
        ax2 = fig.add_subplot(gs[0, 0])
        ax3 = fig.add_subplot(gs[1, 1])
        #ax4 = fig.add_subplot(gs[0, 1])
        #ax5 = fig.add_subplot(gs[0, 2])
        #ax6 = fig.add_subplot(gs[1, 2])

        #rho_gridded = np.copy(aa.rho)
        #rho_gridded[:,:,0,:] = 1000

        #rho_griddedp = np.copy(aa.rho)
        #rho_griddedp[:,0,:,:] = 1000
        #aa.rho[:,:,0,:] = 10000
        #aa.rho[:,0,:,:] = 10000
        #aa.rho[:,:,:,0] = 10000
 
        aa.get_potentials(get_accel=True, get_companion_grav=True)        
        tidal = -1 * aa.rho* aa.gradient(aa.companion_grav_pot + aa.accel_pot, coordinates=grid_type)
        tidal_c = aa.get_C_vec(tidal)
        tidal_c = aa.native_to_cart(tidal_c)
        tidal_c_mag = get_magnitude(tidal_c)

        lrl_native, lrl_cart = aa.get_lrl(components = True)        
        total_mass = aa.integrate(aa.rho, variable='all')
        rhoxlrl = [aa.rho*lrl_cart[0], aa.rho*lrl_cart[1], 0]
        total_rhoxlrl = aa.vec_vol_integrate(rhoxlrl)

        pos_x_lrl_bool = (rhoxlrl[0] >= 0)
        pos_x_lrl_y = np.multiply(pos_x_lrl_bool, rhoxlrl[1])
        neg_x_lrl_bool = np.logical_not(pos_x_lrl_bool)
        neg_x_lrl_y = np.multiply(neg_x_lrl_bool, rhoxlrl[1])
        lrl_orient = (np.arctan(pos_x_lrl_y / rhoxlrl[0]) % (2*np.pi)) + np.multiply(((np.pi + np.arctan(neg_x_lrl_y / rhoxlrl[0])) % (2*np.pi)), neg_x_lrl_bool)

        mass_weighted_eccent = total_rhoxlrl / total_mass

        aa.get_primaries(get_press=True)

        aa.midplane_colorplot(aa.press, ax2, slicetype='z', log=True, vbound=[1e-5,1e2])
        aa.midplane_colorplot(aa.rho, ax1, slicetype='z', log=True, vbound=[1e-5,1e2])
        ax3.plot(aa.possible_r, aa.integrate(aa.press, "shell"))
        #aa.midplane_colorplot(tidal_c[1], ax3, slicetype='z', log=False, vbound=[-1e0,1e0])
        #aa.midplane_colorplot(lrl_orient, ax4, slicetype='z', log=False, vbound=[0,2], angular=True) 
        #aa.midplane_colorplot(aa.rho, ax2, log=True, vbound=[1e-5,1e1])
        #aa.midplane_colorplot(drhodr, ax3, slicetype='y', log=False, vbound=[-2,2])
        #aa.midplane_colorplot(drhodz, ax6, slicetype='y', log=False, vbound=[-2,2])
        #aa.midplane_colorplot(drhodr, ax4, log=False, vbound=[-2,2])
        #aa.midplane_colorplot(drhodphi, ax5, log=False, vbound=[-2,2])

        ax2.set_title(r"Pressure")
        ax1.set_title(r"Density")
        ax3.set_title(r"Pressure")
        #ax6.set_title(r"$\frac{d\rho}{dz}$")
        #ax4.set_title(r"mass rescaled lrl")
        #ax5.set_title(r"$\frac{d\rho}{d\phi}$")
        
        #ax1.set_yticks(ticks=[-12, -6, 0, 6, 12], labels=[-0.5, -0.25, 'test', 0.25, 0.5])

        plt.tight_layout()
        plt.savefig("%s/%s%s_test_plot%05d.png" % (savedir, dname, aname, fnum))
        plt.close()


#main()

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

ff2_test()