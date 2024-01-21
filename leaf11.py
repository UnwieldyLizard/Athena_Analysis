from branches.angular_momentum import *
#from branches.momentum import *
#from branches.B_fields import *
#from branches.meccentricity_plot import*
#from branches.meccentricity_profile import *
#from branches.profile import *
#from branches.roots.misc_func import *
#from branches.roots.utility import *
#from branches.mplottesting import *

def simple_loop(dname, fnum_range, file_spacing, function):
    now = datetime.now()
    fnum_range = np.arange(fnum_range[0], fnum_range[-1]+1, file_spacing)
    for fnum in fnum_range:
        function(dname, fnum)
        try:
            function(dname, fnum)
        except:
            logging.info(f"operation failed on fnum = {fnum}, exiting... (Don't panic this probably just means you've gone through all the data)")
            break

#main("Cyl_14", [0, 10000])
#eccentricity_profile("Cyl_13_Stream", [0,10000])


#simple_loop("Cyl_11", [0, 10000], 500, eccentricity_radial_profile)

#B = B_fields("Cyl_15")
#B.vector_profile(1200)
#B.alfven(1200)

#B = B_fields("Cyl_14")
#B.vector_profile(3000)
#B.alfven(3325)

#B = B_fields("Cyl_11")
#B.alpha_plot(3325)
#B.vector_profile(3325)
#B.profile(0)

#A = AngularMomentum("Cyl_7")
#A.profile(5001, vbound = 1e5)

#A = AngularMomentum("Cyl_7")
#A.profile(4000, flatten=True, just_torques=True)
#del A
#A = AngularMomentum("Cyl_13_Stream")
#A.profile(1800, flatten=True, just_torques=True)
#del A
A = AngularMomentum("Cyl_14")
A.profile(1350, flatten=True, just_torques=True)

#M = Momentum("Cyl_14")
#M.evolution([0,10000])

#eccentricity_plot("Cyl_15", [0, 10000], inertial=True)
#compare_accretion(["Cyl_11", "Cyl_12", "Cyl_13", "Cyl_14", "Cyl_15"], [131, 10000])#, restart=True)
#U = Utility("Cyl_14", "_vertKE_test_plot")
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