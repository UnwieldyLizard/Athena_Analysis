#from branches.angular_momentum import *
#from branches.momentum import *
#from branches.B_fields import *
from branches.meccentricity_profile import*
#from branches.profile import *
#from branches.roots.misc_func import *
#from branches.roots.utility import *
#from branches.mplottesting import *

#main("Cyl_13", [0, 10000])

#M = Momentum("Cyl_13")
#M.evolution([0,10000])
#B = B_fields("Cyl_13")
#B.alpha_beta_radial(2200)

split_profile("Cyl_15_2", [0, 1000000 ])

#eccentricity_plot("Cyl_15", [0, 10000], inertial=True)
#compare_accretion(["Cyl_11", "Cyl_12", "Cyl_13", "Cyl_14", "Cyl_15"], [131, 10000])#, restart=True)

#U = Utility("Cyl_13", "_vertKE_test_plot")
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