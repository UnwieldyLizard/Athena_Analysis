#from branches.angular_momentum import *
from branches.meccentricity_plot import*
#from branches.profile import *
from branches.roots.misc_func import *
#from branches.roots.utility import *

eccentricity_plot("Cyl_12", [0, 10000], inertial=True)
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