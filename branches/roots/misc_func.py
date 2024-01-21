from .athena_analysis import *
import numpy as np
import logging
from datetime import datetime

def mkdir_if_not_exist(path):
    if not os.path.exists(path):
        os.makedirs(path)

def find_file_spacing(dname, start_point = 0):
    data_location = file.data_loc + dname
    max_fnum = 1 + start_point
    file_spacing = None
    while file_spacing is None:
        while not os.path.exists("%s/disk.out1.%05d.athdf" % (data_location, max_fnum)):
            if max_fnum > 50:
                raise("didn't find any of first 50 files, I suspect your file path is wrong")
            max_fnum += 1
        file_spacing = max_fnum
    return max_fnum

def simple_loop(fnum_range, file_spacing, function):
    now = datetime.now()
    fnum_range = np.arange(fnum_range[0], fnum_range[-1]+1, file_spacing)
    for fnum in fnum_range:
        logging.info(datetime.now()-now)
        now = datetime.now()
        logging.info("fnum = %d" % fnum)
        try:
            function(fnum)
        except:
            logging.info(f"operation failed on fnum = {fnum}, exiting... (Don't panic this probably just means you've gone through all the data)")
            break

def dirty_loop(fnum_range, file_spacing, function, dname, aname, sname=""):
    """
    Intended for dealing with partial data sets or filling in gaps of incomplete jobs.
    This loop skips timesteps where it detects and output already exists or data is missing.
    """
    now = datetime.now()
    fnum_range = np.arange(fnum_range[0], fnum_range[-1]+1, file_spacing)
    for fnum in fnum_range:
        logging.info(datetime.now()-now)
        now = datetime.now()
        logging.info("fnum = %d" % fnum)
        if os.path.exists(file.savedir+dname+"/"+dname+aname+"/"+dname+aname+"%05d%s.png" % (fnum, sname)):
            logging.info("Output already exist, skipping")
            continue

        if not os.path.exists(file.data_loc + dname + "/disk.out1.%05d.athdf" % (fnum)):
            logging.info("'" + file.data_loc + dname + "/disk.out1.%05d.athdf' does not exist" % (fnum))
            logging.info("Since data file missing, skipping")
            continue
        
        function(fnum)


def radial_slice_loop(dname, fnum, grid_type, function, scale_factor=1, start_idx=0):
    now = datetime.now()
    dummy_instance = Athena_Analysis("%s/disk.out1.%05d.athdf" % (file.data_loc + dname, fnum), grid_type=grid_type)
    dummy_instance._axes()
    phi_list = np.zeros(int(len(dummy_instance.possible_phi) / scale_factor))
    for p in range(len(phi_list)):
        phi_list[p] = dummy_instance.possible_phi[scale_factor * p]
    del dummy_instance
    
    for i, phi in enumerate(phi_list):
        if start_idx > len(phi_list):
            logging.info("starting index beyond range of angle indices, proceeding to pegging")
            break
        if i < start_idx:
            if i == (start_idx-1):
                logging.info(f"skipped until: i = {start_idx}")
            continue
        logging.info(datetime.now()-now)
        now = datetime.now()
        logging.info(f"phi = {phi:.2f}"+", %d/%d, i=%d" % (i+1, len(phi_list), i))
        function(dname, fnum, grid_type, phi, radial_slice_loop=True)
        #renaming
        if i == 0 or i == start_idx:
            old_name_base = "%s/%s%s%05d_phi=" % (file.savedir + dname + "/" + dname + function.aname, dname, function.aname, fnum) 
            new_name_base = "%s/%s%s%05d" % (file.savedir + dname + "/" + dname + function.aname, dname, function.aname, fnum)
        old_name = old_name_base+f"{phi/np.pi:.2f}pi"+".png"
        new_name = new_name_base+"_angular_idx:%05d.png" % i
        os.rename(old_name, new_name)

    #pegging
    framerate = int(len(phi_list) / 60)
    if framerate == 0:
        logging.info("Lmao this'll be a shitty video with so few frames, setting framerate to 1.")
        framerate = 1
    os.system("ffmpeg -framerate "+str(framerate)+" -i "+new_name_base+"_angular_idx:%05d.png -vf scale=1280:-2 "+new_name_base+"_angular_loop"+".mp4")
