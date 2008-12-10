#!/usr/bin/env python
import os
import shutil
import subprocess
import time
import string

Ens=2



WINDS="/mnt/ocg0/DISK2/abarth/WFSROMS/Data/Ens_2005/"
WINDS="/mnt/ocgmod3/home/abarth/WFSROMS/Data/OIWinds/Ens_2005/"

model="/home/abarth/WFSROMS/WFS/roms-2.0_nr/oceanS"

f=open('wfs.assim.codar_wind.INIT','r')

for i in f:
  if i.find('ErrorSpace.dimension') != -1:
    Ens=string.atoi(i.split('=')[1])

f.close()

print "Ens ",Ens

def removefile(file):
  try:
    os.remove(file)
  except:
    pass
      
    

for j in range(1,Ens+1):  
  dir="%03d" % (j)
  if not os.path.isdir(dir):
    os.mkdir(dir)
  shutil.copyfile("WFS_ini.001.nc",dir + "/WFS_ini.nc");
  removefile(dir + "/WFS.out")


removefile("WFS.out")
removefile("assim.log")

process = []

max_running=8

for i in range(0,50):

  members = range(Ens,0,-1)
  process = [];
  
  while len(members) != 0 or len(process) != 0:
      while len(process) < max_running and len(members) > 0:
          # launch new
          j = members.pop()-1
          print "member ",j+1
          dir="%03d" % (j+1)
          os.chdir(dir)
          removefile("wind_oi_ens.nc")
          os.symlink(WINDS + "/wind_oi_ens%03d.nc" % (j+1),"wind_oi_ens.nc")

          cmd=model + " < ../Init/WFS.%03d.in >> WFS.out" % (i+1)
          print cmd, "in ",os.getcwd()   
          process.append(subprocess.Popen(cmd,shell=True));

          os.chdir("..")

      # look for finished process

      for p in process:
        if p.poll() != None:
#          print "finished"
          process.remove(p)

      time.sleep(5)

  # Launch the assimilation
  
  subprocess.Popen("./assim_ens.sh wfs.assim.codar_wind.INIT %d" % (i+1),shell=True).wait()



     
     

  
    
    
    
