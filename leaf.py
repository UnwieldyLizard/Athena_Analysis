#from branches.angular_momentum import *
#from branches.comparison import *
from branches.meccentricity_profile import *
#from branches.eccentricity import*
#from branches.profile import *
#from branches.roots.misc_func import *
#from branches.roots.utility import *
#from branches.momentum import *

#eccentricity_radial_profile("Cyl_11", 2000)

"""
def simple_loop(fnum_range, file_spacing, function):
    now = datetime.now()
    fnum_range = np.arange(fnum_range[0], fnum_range[-1]+1, file_spacing)
    for fnum in fnum_range:
        try:
            function(fnum, flatten=True, just_torques=True)
        except:
            logging.info(f"operation failed on fnum = {fnum}, exiting... (Don't panic this probably just means you've gone through all the data)")
            break
"""

#A = AngularMomentum("Cyl_13_Stream")
#simple_loop([1650, 1950], 1, A.profile)

#U = Utility("Cyl_13", "_angular_momentum", "spread")
#U.tar()
#del U

#A = AngularMomentum("Cyl_13")

#simple_loop([1500, 2000], 100, A.profile)

#C = Comparison(["Cyl_11", "Cyl_12", "Cyl_13", "Cyl_13_Stream", "Cyl_14", "Cyl_15"])
#C.alpha([0, 4000])

#compare_accretion(["Cyl_11", "Cyl_12", "Cyl_13", "Cyl_13_Stream", "Cyl_14", "Cyl_15"], [0, 100000], restart = True, plot_every = 100, pickle_every = 100)
#append_dataset_accretion(["Cyl_13_Stream"], [0,10000])

#momentum_plot("Cyl_11", [0, 10000], bar_widths=[20,30,5])

#compare_alpha(["Cyl_11", "Cyl_12", "Cyl_13", "Cyl_14", "Cyl_15"], [1600, 10000], res=True, B_field_focus=True)

#eccentricity_profile("Cyl_13_2", [0, 100000])
#replot("Cyl_11_2", 3400)
replot("Cyl_15_2", 3100, "_recomp", stress_test=True, recompute_sum=True)

#E = Eccentricity("Cyl_11_2")
#dirty_loop([0, 3000], 1, E.plot, "Cyl_11_2", "_eccent")

#U = Utility("Cyl_11_2", "_eccent")
#U.tar()
#eccentricity_plot("Cyl_15", [409, 10000], inertial=True)
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
