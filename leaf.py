#from branches.angular_momentum import *
#from branches.meccentricity_plot import*
from branches.profile import *
from branches.roots.misc_func import *
from branches.roots.utility import *

compare_beta(["Cyl_11", "Cyl_12", "Cyl_13"], [556, 10000])

"""
P = Profile("Cyl_11")
P.plasma_beta(225)
del P
P = Profile("Cyl_12")
P.plasma_beta(225)
del P
P = Profile("Cyl_13")
P.plasma_beta(225)
del P
"""