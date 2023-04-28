from branches.angular_momentum import *
from branches.roots.misc_func import *
from branches.roots.utility import *

L = AngularMomentum("Cyl_1")
L.profile(1000)
L.profile(2000)
L.profile(3000)
L.profile(4000)
L.profile(5000)
#L.time_average_profile(500, 1)