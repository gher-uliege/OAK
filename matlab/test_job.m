
job = Job({'/u/abarth/bin/run_octave.sh ls'});
isrunning(job)
wait(job)
disp('finished')
