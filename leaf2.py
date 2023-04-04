#from roots.misc_func import *
#from meccentricity_profile import *
#tidal_profile("Cyl_1", 5000, "Cylindrical")
#radial_slice_loop("Cyl_1", 2000, "Cylindrical", tidal_profile, start_idx=506)
#del reduce_dimensions.normalization_weight
#simple_loop("Cyl_1", [2090,10000], 10, "Cylindrical", tidal_profile)

from branches.waves import *
#fourier_waves_loop("Cyl_1", [5200, 10000], 10, [0,5],"Cylindrical")
dname = "Cyl_7"
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