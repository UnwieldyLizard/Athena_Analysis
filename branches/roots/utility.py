import os
import logging
from .file_config import*
file = File()

logging.basicConfig(filename=file.logs_loc+"/utility.log", encoding='utf-8', level=logging.INFO)

class Utility():
    def __init__(self, dname, aname):
        self.dname = dname
        self.aname = aname

    def renumber(self):
        for i, n in enumerate(range(3093, 10000, 4)):
            old_name = file.savedir+self.dname+"/"+self.dname+self.aname+"/"+self.dname+self.aname+"%05d.png" % n
            new_name = file.savedir+self.dname+"/"+self.dname+self.aname+"/"+self.dname+self.aname+"%05d.png" % (i+773)
            if os.path.exists(old_name):
                os.rename(old_name, new_name)
            else:
                print(f"renumbered up to n={n-1}, i={i-1}")
                break

    def peg(self):
        os.system("ffmpeg -framerate 12 -i "+file.savedir+self.dname+"/"+self.dname+self.aname+"/"+self.dname+self.aname+"%05d.png -vf scale=1280:-2 "+file.savedir+self.dname+"/"+self.dname+self.aname+"/"+self.dname+self.aname+".mp4")

    def tar(self):
        os.system("tar -czvf "+file.savedir+self.dname+"/"+self.dname+self.aname+"/"+self.dname+self.aname+".tar.gz "+file.savedir+self.dname+"/"+self.dname+self.aname)#+"/"+self.dname+self.aname)

    def radial_peg(self, fnum):
        new_name_base = "%s/%s%s%05d" % (file.savedir + self.dname + "/" + self.dname + self.aname, self.dname, self.aname, fnum)
        max_i = 0
        while os.path.exists(new_name_base+f"_angular_idx_{max_i:05}.png"):
            max_i += 1
        logging.info(max_i)
        framerate = int(max_i / 60)
        if framerate == 0:
            logging.info("Lmao this'll be a shitty video with so few frames, setting framerate to 1.")
            framerate = 1
        os.system("ffmpeg -framerate "+str(framerate)+" -i "+new_name_base+"_angular_idx_%05d.png -vf scale=1280:-2 "+new_name_base+"_angular_loop"+".mp4")