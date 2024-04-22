#from branches.angular_momentum import *
from branches.waves import *
#from branches.momentum import *
#from branches.meccentricity_plot import *
#from branches.meccentricity_profile import *
#from branches.profile import *
#from branches.roots.misc_func import *
#from branches.roots.utility import *
#from branches.B_fields import *

F = Fourier_Waves("Cyl_11_2")
#F.collect_density(start_orbit=60)
#F.collect_momentum_density(start_orbit=60)
#F.collect_rholrl(start_orbit=60)
#F.plot_momentum_density(fnum_range=[2100, 2205], r_slicepoint=sim.three_one_res, bounds=[0, 1.5])
F.tidal_waves(fnum_range=[2100, 2205], wave_modes=[3,2,1,0], sname="")

#momentum_plot("Cyl_14", [0, 10000], bar_widths=[20,30,5])

#B = B_fields("Cyl_13")
#simple_loop([0, 10000], 1, B.profile)

#B = B_fields("Cyl_14")
#B.evolution([0, 10000])

#eccentricity_plot("Cyl_14", [0, 10000], inertial=True)
#split_profile_rates("Cyl_11", 1600, "rates", num_orbits_ave=1)
#split_profile_rates("Cyl_12", 2299, "rates", num_orbits_ave=1)
#split_profile_rates("Cyl_13", 2200, "rates", num_orbits_ave=1)
#get_late_betas(["Cyl_11", "Cyl_12", "Cyl_13"], cutoff_orbit=35, load_point=3320)
#U = Utility("Cyl_11", "_eccent")
#U.tar()
#del U
#U = Utility("Cyl_14", "_momentum")
#U.tar()
#eccentricity_plot("Cyl_14", [0, 10000], inertial=True)
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