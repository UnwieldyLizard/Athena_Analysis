#from branches.angular_momentum import *
#from branches.momentum import *
#from branches.B_fields import *
#from branches.meccentricity_profile import *
from branches.thermal import *
#from branches.eccentricity import *
#from branches.roots.misc_func import *
#from branches.roots.utility import *
#from branches.profile import *
#from branches.mplottesting import *

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

#A = AngularMomentum("Cyl_13")
#simple_loop([1650, 1950], 1, A.profile)

T = Thermal("Cyl_11_2")
T.temperature_profile(3000)

#replot("Cyl_11_2", 3400, "_recomp", stress_test=True, recompute_sum=True)

#E = Eccentricity("Cyl_11_2")
#dirty_loop([0, 10000], 1, E.plot, "Cyl_15_2", "_eccent")

#B = B_fields("Cyl_11_2")
#B.evolution([0, 10000])
#B.profile(3450)
#B.vector_profile(3450)
#simple_loop([218, 10000], 1, B.beta)

#eccentricity_plot("Cyl_12", [0, 10000], inertial=True)
#U = Utility("Cyl_15_2", "_eccent")
#U.tar()
#split_profile_rates("Cyl_11", 3300, "rates")
#compare_alpha(["Cyl_11", "Cyl_12", "Cyl_13", "Cyl_14"], [0, 10000], res=True)
#append_dataset_alpha(["Cyl_14", "Cyl_15"], [0,10000], res=True, pickle_every=10)
