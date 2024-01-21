#from branches.angular_momentum import *
#from branches.momentum import *
from branches.thermal import *
#from branches.B_fields import *
#from branches.meccentricity_plot import*
#from branches.profile import *
#from branches.roots.misc_func import *
from branches.roots.utility import *

#momentum_plot("Cyl_15", [0, 10000], bar_widths=[20,30,5])

#B = B_fields("Cyl_13_Stream")
#simple_loop([0, 10000], 1, B.pressure_profile)

#B = B_fields("Cyl_15")
#B.evolution([0, 10000])

T = Thermal("Cyl_14")
T.entropy_profile(1950)
T.entropy_profile(2000)

#eccentricity_plot("Cyl_15", [0, 10000], inertial=True)
#compare_accretion(["Cyl_11", "Cyl_12", "Cyl_13", "Cyl_14", "Cyl_15"], [131, 10000])#, restart=True)

#U = Utility("Cyl_13_Stream", "_B_fields", "pressure")
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