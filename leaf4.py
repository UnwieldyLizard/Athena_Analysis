#from branches.comparison import *
#from branches.angular_momentum import *
#from branches.momentum import *
from branches.B_fields import *
#from branches.meccentricity_profile import *
#from branches.eccentricity import *
#from branches.profile import *
#from branches.roots.misc_func import *
#from branches.roots.utility import *
#from branches.mplottesting import *

B = B_fields("Cyl_11_2")
B.vector_profile(3000)
#simple_loop([0, 10000], 1, B.pressure_profile)

#A = AngularMomentum("Cyl_7")
#A.profile(4000, flatten = True, just_torques=True, spread=False)
#A.profile(4000, flatten = True, just_torques=True)

#simple_loop([0, 10000], 1, A.profile)

#momentum_plot("Cyl_13", [0, 10000], bar_widths=[20,30,5])

#B = B_fields("Cyl_12")
#B.evolution([0, 10000])

#E = Eccentricity("Cyl_13_aB", 2201)
#dirty_loop([2201, 3000], 1, E.plot, "Cyl_13_aB", "_eccent")

#eccentricity_profile("Cyl_11_2", [0, 100000], 10, stress_test = True)

#U = Utility("Cyl_13_Stream", "_B_fields", "pressure")
#U.tar()

#split_profile_rates("Cyl_12", 2999, "rates")
#append_dataset_alpha(["Cyl_14", "Cyl_15"], [0,10000])
#replot_compare_alpha(2200, res=True)

#eccentricity_plot("Cyl_13", [0, 10000], inertial=True)
#main("Cyl_13", [300, 301])#, 100)
#main("Cyl_13", [703, 704])#, 100)

#P = Profile("Cyl_11_2")
#P.temporal_profile([0, 100000])

#C = Comparison(["Cyl_7", "Cyl_11_2", "Cyl_13_2", "Cyl_15_2"], 4)
#C.beta([0,5000])
#C.alpha_replot("full", log=False, ylims=[[-0.01,0.15],[-5,15],[-0.2,1]])# ylims=[0,15])