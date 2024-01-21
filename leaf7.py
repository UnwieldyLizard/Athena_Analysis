#from branches.angular_momentum import *
from branches.momentum import *
#from branches.energy import *
#from branches.meccentricity_plot import*
#from branches.profile import *
#from branches.roots.misc_func import *
#from branches.roots.utility import *
#from branches.mplottesting import *

M = Momentum("Cyl_11")
M.evolution([0,10000])

#kinetic_energy_plot(["Cyl_11", "Cyl_12", "Cyl_13", "Cyl_13_Stream", "Cyl_14", "Cyl_15"], [0, 10000])

#main("Cyl_11", [0, 10000])

#eccentricity_plot("Cyl_12", [2327, 10000], inertial=True)
#U = Utility("Cyl_11", "_vertKE_test_plot")
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