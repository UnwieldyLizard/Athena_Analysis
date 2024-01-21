#from branches.angular_momentum import *
from branches.momentum import *
#from branches.meccentricity_plot import*
#from branches.profile import *
#from branches.roots.misc_func import *
from branches.roots.utility import *
#from branches.mplottesting import *

#main("Cyl_12", [0, 10000])

M = Momentum("Cyl_12")
M.evolution([0,10000])

# from branches.roots.athena_analysis import *
# data_location = file.data_loc + "Cyl_13"
# filename = "%s/disk.out1.%05d.athdf" % (data_location, 0)
# aa = Athena_Analysis(filename=filename, grid_type=file.grid_types["Cyl_13"])
# aa.get_Bfields()
# b_z = aa.integrate(aa.B_z, "phi", "Full", "r")
# logging.info("Unique B_z Values")
# logging.info(set((aa.B_z).flatten()))
# logging.info("B_z integrated through r and phi")
# logging.info(b_z)

#eccentricity_plot("Cyl_15", [0, 10000], inertial=True)
#compare_accretion(["Cyl_11", "Cyl_12", "Cyl_13", "Cyl_14", "Cyl_15"], [131, 10000])#, restart=True)
#U = Utility("Cyl_12", "_vertKE_test_plot")
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