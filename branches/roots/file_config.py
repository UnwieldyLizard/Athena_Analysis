from pathlib import Path

class File:
    def __init__(self):
        #save point and data sources
        self.athena_analysis_directory = str(Path(__file__).parents[2])
        self.savedir = "/mnt/c/Users/morga/Desktop/School/research_stuff/processed_data/"
        self.data_loc_D = "/mnt/d/Research_Data/"
        self.data_loc = "/mnt/e/Research_Data/"
        self.logs_loc = self.athena_analysis_directory + "/branches/logs"
        self.mkitp = "/home/morgan/mnt/kitp/"
        self.medd = "/home/morgan/mnt/edd/"
        #The data set bryance left us
        self.cvthin2 = "data2/cvdisk/CVThin2/Data"
        self.cvthick3 = "data/cvdisk/CVThick3"
        self.alpha3 = "data/cvdisk/superhump_3d_alpha03"

        #set up dictionaries
        self.grid_types = {}
        self.grid_types["Cyl_1"] = "Cylindrical"
        self.grid_types["Cyl_2"] = "Cylindrical"
        self.grid_types["Cyl_6"] = "Cylindrical"
        self.grid_types["Cyl_7"] = "Cylindrical"

5