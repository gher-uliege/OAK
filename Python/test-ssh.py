import time
import subprocess

workdir='/home/abarth/Test/'
model="sleep 10"


def run_octave(cmd):
    p = subprocess.Popen(['octave','-q'],stdin=subprocess.PIPE)
    p.stdin.write(cmd)
    # stops octave
    p.stdin.close()
    p.wait()
    

cmd = """
A = [1 2 3 4]
"""

run_octave(cmd)
