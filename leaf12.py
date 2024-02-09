#from branches.angular_momentum import *
#from branches.B_fields import *
#from branches.eccentricity import *
from branches.meccentricity_profile import *
#from branches.profile import *
#from branches.roots.misc_func import *
#from branches.roots.utility import *
#from branches.mplottesting import *

#main("Cyl_15", [0, 10000])

eccentricity_profile("Cyl_15_2", [3101, 100000], stress_test=True)

#M = Momentum("Cyl_15")
#M.evolution([0,10000])

#B = B_fields("Cyl_14")
#B.vector_profile(1)
#B.evolution([0,10000])
#B.vector_profile(3)

#E = Eccentricity("Cyl_14")
#dirty_loop([0, 10000], 1, E.plot, "Cyl_14", "_eccent")

#eccentricity_plot("Cyl_15", [0, 10000], inertial=True)
#compare_accretion(["Cyl_11", "Cyl_12", "Cyl_13", "Cyl_14", "Cyl_15"], [131, 10000])#, restart=True)
#U = Utility("Cyl_14", "_eccent")
#U.tar()
"""
P = Profile("Cyl_11")
P.plasma_beta(1500)
del P
P = Profile("Cyl_12")
P.plasma_beta(1500)
del P
P = Profile("Cyl_13")
P.plasma_beta(1500)
del P
"""