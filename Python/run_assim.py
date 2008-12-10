#!/usr/bin/env python
import os
import shutil
import subprocess
import time
import string

WINDS="/mnt/ocg0/DISK2/abarth/WFSROMS/Data/Ens_2005/"
WINDS="/mnt/ocgmod3/home/abarth/WFSROMS/Data/OIWinds/Ens_2005/"



nodes = [{"hostname": "ocgmod1.marine.usf.edu",
          "executable": "/mnt/ocgcluster/home/abarth/WFSROMS/WFS/roms-2.0_nr/oceanS",
          "num_threads": 1,
          "workdir": '/mnt/ocgcluster/home/abarth/WFSROMS/Output_assim_codar_filter_s12_evolutive/',
          },         
         {"hostname": "ocgmod3.marine.usf.edu",
          "executable": "/mnt/ocgcluster/home/abarth/WFSROMS/WFS/roms-2.0_nr/oceanS",
          "num_threads": 4,
          "workdir": '/mnt/ocgcluster/home/abarth/WFSROMS/Output_assim_codar_filter_s12_evolutive/',
          },
         {"hostname": "ocg9.marine.usf.edu",
          "executable": "/mnt/ocgcluster/home/abarth/WFSROMS/WFS/roms-2.0_nr/oceanS",
          "num_threads": 4,
          "workdir": '/mnt/ocgcluster/home/abarth/WFSROMS/Output_assim_codar_filter_s12_evolutive/',
          },
         {"hostname": "ocgcluster.marine.usf.edu",
          "executable": "/mnt/ocgcluster/home/abarth/WFSROMS/WFS/roms-2.0_nr/oceanS",
          "num_threads": 8,
          "workdir": '/mnt/ocgcluster/home/abarth/WFSROMS/Output_assim_codar_filter_s12_evolutive/',
          }]

##nodes = [{"hostname": "ocgmod3.marine.usf.edu",
##          "executable": "/mnt/ocgcluster/home/abarth/WFSROMS/WFS/roms-2.0_nr/oceanS",          
##          "num_threads": 4,
##          "workdir": '/mnt/ocgcluster/home/abarth/WFSROMS/Output_assim_codar_filter_s12_evolutive/',
##          }]

#          "input": "../Init_debug/WFS.%03d.in",

config = {
          "input": "../Init/WFS.%03d.in",
          "output": "WFS.%03d.out",
          "workdir": '/mnt/ocgcluster/home/abarth/WFSROMS/Output_assim_codar_filter_s12_evolutive/',
}

def removefile(file):
  try:
    os.remove(file)
  except:
    pass

def run_ensemble(nodes,config,Ens,i):
  members = range(Ens,0,-1)

  for n in nodes:
    n['process'] = [];
  
  while len(members) != 0 or running != 0:
    for n in nodes:    
      while len(n['process']) < n['num_threads'] and len(members) > 0:
          # launch new
          j = members.pop()-1
          #print "member ",j+1
          directory=n['workdir'] + "%03d" % (j+1)
          inputfile=config['input'] % (i+1)
          outputfile=config['output'] % (i+1)
          
          cmd = 'cd %s; %s < %s > %s' % (directory,n['executable'],inputfile,outputfile)

          print "running ",cmd," on ",n['hostname']

          p = subprocess.Popen(['ssh',n['hostname'],cmd])
          n['process'].append(p);

      # look for finished process

    running = 0
      
    for n in nodes:    
      for p in n['process']:
        if p.poll() != None:
          print "finished: return code: ",p.poll()
          #print p.stdout.read()
          n['process'].remove(p)
          
      running = running + len(n['process'])

    #print 'running ',running          

    time.sleep(1)

def prep_envrionment(Ens):
  for j in range(1,Ens+1):  
    dir="%03d" % (j)
    if not os.path.isdir(dir):
      os.mkdir(dir)
    #shutil.copyfile("WFS_ini.001.nc",dir + "/WFS_ini.nc");
    shutil.copyfile("/mnt/ocg9/home/abarth/WFSROMS/Output_wfs_ens_100/Ens%03d/WFS_ini.nc" % j,dir + "/WFS_ini.nc");
    removefile(dir + "/WFS.out")
    removefile(dir + "/wind_oi_ens.nc")
    os.symlink(WINDS + "/wind_oi_ens%03d.nc" % (j+1),dir + "/wind_oi_ens.nc")
    
  removefile("WFS.out")
  removefile("assim.log")

    

def run_assim(i):
  print "assimilation ",i
  subprocess.Popen("./assim_ens.sh wfs.assim.codar_wind.INIT %d" % (i+1),shell=True).wait()


def getInitValue(filename,key):
  f=open(filename,'r')

  for i in f:
    if i.find(key,) != -1:
      val=string.atoi(i.split('=')[1])

  f.close()

  return val

def run_octave(cmd):
    print cmd
    p = subprocess.Popen(['octave','-q'],stdin=subprocess.PIPE)
    p.stdin.write(cmd)
    # stops octave
    p.stdin.close()
    p.wait()

def post_ensemble(Ens):
    run_octave("post_ens(%d)" % Ens);

def post_assim(Ens,i):
    run_octave("post_assim(%d,%d)" % (Ens,i+1));


  
if __name__ == '__main__':
  Ens = getInitValue('wfs.assim.codar_wind.INIT','ErrorSpace.dimension')
  print "Ens",Ens

  for n in nodes:
    n['workdir'] = '/mnt/ocgcluster/' + os.getcwd() + '/'
  
  prep_envrionment(Ens)

  # loop over time
  for i in range(0,50):
    run_ensemble(nodes,config,Ens,i)
    post_ensemble(Ens)
    run_assim(i)
    post_assim(Ens,i)




     
     

  
    
    
    
