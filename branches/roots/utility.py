import os
import logging
from .file_config import*
file = File()

logging.basicConfig(filename=file.logs_loc+"/utility.log", encoding='utf-8', level=logging.INFO)

class Utility():
    def __init__(self, dname, aname, sname=""):
        self.dname = dname
        self.aname = aname
        self.sname = sname
        self.file_spacing = None
        max_fnum = 1
        while self.file_spacing is None:
            while not os.path.exists(file.savedir+self.dname+"/"+self.dname+self.aname+"/"+self.dname+self.aname+"%05d%s.png" % (max_fnum, self.sname)):
                max_fnum += 1
            self.file_spacing = max_fnum

    def renumber(self):
        for i, n in enumerate(range(0, 10000, self.file_spacing)):
            old_name = file.savedir+self.dname+"/"+self.dname+self.aname+"/"+self.dname+self.aname+"%05d%s.png" % (n, self.sname)
            new_name = file.savedir+self.dname+"/"+self.dname+self.aname+"/"+self.dname+self.aname+"%05d%s.png" % (i, self.sname)
            if os.path.exists(old_name):
                os.rename(old_name, new_name)
            else:
                logging.info(f"renumbered up to n={n-1}, i={i-1}")
                self.file_spacing = 1
                break

    def peg(self):
        if self.file_spacing == 1:
            os.system("ffmpeg -framerate 12 -i "+file.savedir+self.dname+"/"+self.dname+self.aname+"/"+self.dname+self.aname+"%05d"+self.sname+".png -vf scale=1280:-2 "+file.savedir+self.dname+"/"+self.dname+self.aname+"/"+self.dname+self.aname+self.sname+".mp4")
        elif self.file_spacing == 10:
            os.system("ffmpeg -framerate 12 -i "+file.savedir+self.dname+"/"+self.dname+self.aname+"/"+self.dname+self.aname+"%04d0"+self.sname+".png -vf scale=1280:-2 "+file.savedir+self.dname+"/"+self.dname+self.aname+"/"+self.dname+self.aname+self.sname+".mp4")

    def tar(self):
        tar_name = file.savedir+self.dname+"/"+self.dname+self.aname+"/"+self.dname+self.aname+".tar.gz"
        if os.path.exists(tar_name):
            os.system("rm "+tar_name) 
        os.system("tar -czvf "+tar_name+" "+file.savedir+self.dname+"/"+self.dname+self.aname+self.sname)#+"/"+self.dname+self.aname)

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

    def make_movie(self):
        if (self.file_spacing != 1) and (self.file_spacing != 10):
            self.renumber()
        self.peg()