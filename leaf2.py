from roots.misc_func import *
from meccentricity_profile import *
#tidal_profile("Cyl_1", 5000, "Cylindrical")
#radial_slice_loop("Cyl_1", 2000, "Cylindrical", tidal_profile, start_idx=506)
#del reduce_dimensions.normalization_weight
simple_loop("Cyl_1", [2090,10000], 10, "Cylindrical", tidal_profile)