import os
import logging
from roots.file_config import*
file = File()

logging.basicConfig(filename=file.logs_loc+"/utility.log", encoding='utf-8', level=logging.INFO)

dname = "Cyl_7"
#aname = "_eccent"
aname = "_tidal"

def renumber():
    for i, n in enumerate(range(3093, 10000, 4)):
        old_name = file.savedir+dname+"/"+dname+aname+"/"+dname+aname+"%05d_az_ave.png" % n
        new_name = file.savedir+dname+"/"+dname+aname+"/"+dname+aname+"%05d_az_ave.png" % (i+773)
        if os.path.exists(old_name):
            os.rename(old_name, new_name)
        else:
            print(f"renumbered up to n={n-1}, i={i-1}")
            break

def peg():
    os.system("ffmpeg -framerate 12 -i "+file.savedir+dname+"/"+dname+aname+"/"+dname+aname+"%05d_az_ave.png -vf scale=1280:-2 "+file.savedir+dname+"/"+dname+aname+"/"+dname+aname+"_az_ave.mp4")

def radial_peg(fnum):
    new_name_base = "%s/%s%s%05d" % (file.savedir + dname + "/" + dname + aname, dname, aname, fnum)
    max_i = 0
    while os.path.exists(new_name_base+f"_angular_idx_{max_i:05}.png"):
        max_i += 1
    logging.info(max_i)
    framerate = int(max_i / 60)
    if framerate == 0:
        logging.info("Lmao this'll be a shitty video with so few frames, setting framerate to 1.")
        framerate = 1
    os.system("ffmpeg -framerate "+str(framerate)+" -i "+new_name_base+"_angular_idx_%05d.png -vf scale=1280:-2 "+new_name_base+"_angular_loop"+".mp4")

renumber()
peg()
#radial_peg(2000)