#from branches.angular_momentum import *
#from branches.comparison import *
#from branches.B_fields import *
#from branches.meccentricity_profile import *
#from branches.eccentricity import*
#from branches.profile import *
#from branches.roots.misc_func import *
from branches.roots.utility import *
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

#C = Comparison(["Cyl_7", "Cyl_11_2", "Cyl_13_2", "Cyl_15_2"], 4)
#C.beta([0,5000])
#C.alpha_replot("paper", log=False, ylims=[[-0.01,0.15],[-5,15],[-0.2,1]])

#E = Eccentricity("Cyl_13_2")
#E.plot(3451, sname="test:log")
#simple_loop([0, 10000], 1, E.plot)

U = Utility("Cyl_13_2", "_eccent")
U.tar()

'''
filename = "%s/disk.out1.%05d.athdf" % (file.data_loc + "Cyl_11_2", 0)
aa = Athena_Analysis(filename=filename, grid_type="Cylindrical")
aa.native_grid(get_r=True)
aa.get_Bfields()
aa.get_primaries(get_rho=True)
B = aa.integrate(aa.rho, "shell")
annulus_r = np.ma.masked_where(B < 1e-02, aa.possible_r)
print(aa.possible_rf)
print(annulus_r)
'''

#B = B_fields("Cyl_7")
#B.profile(2000)
#B = B_fields("Cyl_13_2")
#B.alpha_plot(3124)
#B = B_fields("Cyl_15_2")
#B.alpha_plot(4176)

#simple_loop([1500, 2000], 100, A.profile)

#C = Comparison(["Cyl_11", "Cyl_12", "Cyl_13", "Cyl_13_Stream", "Cyl_14", "Cyl_15"])
#C.alpha([0, 4000])

#compare_accretion(["Cyl_11", "Cyl_12", "Cyl_13", "Cyl_13_Stream", "Cyl_14", "Cyl_15"], [0, 100000], restart = True, plot_every = 100, pickle_every = 100)
#append_dataset_accretion(["Cyl_13_Stream"], [0,10000])

#momentum_plot("Cyl_11", [0, 10000], bar_widths=[20,30,5])

#compare_alpha(["Cyl_11", "Cyl_12", "Cyl_13", "Cyl_14", "Cyl_15"], [1600, 10000], res=True, B_field_focus=True)

#eccentricity_profile("Cyl_13_2", [0, 100000])
#replot("Cyl_11_2", 3400)
#replot("Cyl_15_2", 3100, "_recomp", stress_test=True, recompute_sum=True)

#E = Eccentricity("Cyl_11_2")
#dirty_loop([0, 3000], 1, E.plot, "Cyl_11_2", "_eccent")

#U = Utility("Cyl_11_2", "_eccent")
#U.peg()

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
