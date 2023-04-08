from branches.profile import *
from branches.roots.misc_func import *


P = Profile("Cyl_7")
simple_loop([40,5013], 4, P.profile)
#P.profile(fnum=5013, r_slicepoint=2)
#P.profile(fnum=5013, r_slicepoint=4)
#P.profile(fnum=5013, r_slicepoint=6)
#P.profile(fnum=5013, r_slicepoint=8)
#P.profile(fnum=5013, r_slicepoint=10)
#P.profile(fnum=5013, r_slicepoint=12)