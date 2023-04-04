#from roots.misc_func import *
#from meccentricity_profile import *
#replot("Cyl_1", 5215, pname="a_tall_", aspect_ratio=0.5)
#main(dname="Cyl_1", fnum_limits=[5011, 5300], file_spacing=1, plot_every=10, MHD = False, coordinates = "Cylindrical", alpha=0.1)#, restart=True)
#tidal_profile_cmap("Cyl_1", 2000, "Cylindrical")
#simple_loop("Cyl_1", [2880,10000], 10, "Cylindrical", tidal_profile)

from branches.waves import *
#fourier_waves_loop("Cyl_1", [5200, 10000], 10, [0,5],"Cylindrical")
dname = "Cyl_2"
F = Fourier_Waves(dname)
F.collect_data(start_fnum=952)
F.plot()
del F
F = Fourier_Waves(dname)
F.collect_data(start_fnum=1952)
F.plot()
del F
F = Fourier_Waves(dname)
F.collect_data(start_fnum=2952)
F.plot()
del F
F = Fourier_Waves(dname)
F.collect_data(start_fnum=3952)
F.plot()
del F
F = Fourier_Waves(dname)
F.collect_data(start_fnum=4952)
F.plot()

#from orbital import *
#orbital_velocity_analysis("Cyl_1", 5000, "Cylindrical", alpha = 0.1)#, phi_slicepoint=((1)*np.pi/2))
#res_orbit_plot("Cyl_1", 5000, "Cylindrical")
#simple_loop("Cyl_1", [0,10000], 10, "Cylindrical", res_orbit_plot)

#from profile import *
#mass_profile_loop("Cyl_1", [0,10000], 10 , 6, "Cylindrical")