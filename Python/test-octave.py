import time
import subprocess

workdir='/home/abarth/Test/'
model="sleep 10"

cmd = 'cd %s; %s; ls -l' % (workdir,model)

print cmd
p = subprocess.Popen(['octave','--no-history'])
                     stdin=subprocess.PIPE,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE)

p.stdin.write('ls\n')

#p.stdin.write('ls\n')
#p.stdin.write('sleep 10\n')
#p.stdin.write('ls -l\n')
#p.stdin.write('exit\n')
#p.stdin.flush()
#p.stdout.flush()

print p.stdout.readline()
#

while p.poll() == None:
    print "running "
    time.sleep(1)
##    
##
##print p.stdout.read()
