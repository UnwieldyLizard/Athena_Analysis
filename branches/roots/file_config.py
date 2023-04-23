from pathlib import Path

class File:
    def __init__(self):
        #Self and Log locations (YOU DON'T NEED TO TOUCH THIS)
        self.athena_analysis_directory = str(Path(__file__).parents[2])
        self.logs_loc = self.athena_analysis_directory + "/branches/logs"

        #save point and data sources (YOU FILL THIS OUT WITH THE APPROPRIATE FILE PATHS)
        #save_dir should be a general directory for all the results, Athena_Analysis will build and organize all the subdirectories itself.
        self.savedir = "/home/mohana/mhd_superhump/processed_data/"
        #Athena_Analysis will expect the data for each dataset to be in a directory called data_loc+dname
        self.data_loc = "/home/mohana/mhd_superhump/data_raw/globus_dump/"

        #set up dictionaries (YOU FILL THIS OUT)
        #for each item of the dictionary the key should be the dname and the value the grid type
        #grid types should be either "Spherical" or "Cylindrical"
        self.grid_types = {}
        self.grid_types["example_dname"] = "Cylindrical"
        self.grid_types["Cyl_1"] = "Cylindrical"
        self.grid_types["Cyl_2"] = "Cylindrical"
        self.grid_types["Cyl_6"] = "Cylindrical"
        self.grid_types["Cyl_7"] = "Cylindrical"

        self.MHD = {}
        self.MHD["Cyl_1"] = False
        self.MHD["Cyl_2"] = True
        self.MHD["Cyl_7"] = True

        self.alpha = {}
        self.alpha["Cyl_1"] = 0.1
        self.alpha["Cyl_2"] = None
        self.alpha["Cyl_7"] = None

        #my weird shit (I told myself I'd remove this before sharing the github with anyone, if somehow you see this you can safely delete this block of code.)
        self.data_loc_D = "/mnt/d/Research_Data/"
        self.mkitp = "/home/morgan/mnt/kitp/"
        self.medd = "/home/morgan/mnt/edd/"
        #The data set bryance left us
        self.cvthin2 = "data2/cvdisk/CVThin2/Data"
        self.cvthick3 = "data/cvdisk/CVThick3"
        self.alpha3 = "data/cvdisk/superhump_3d_alpha03"