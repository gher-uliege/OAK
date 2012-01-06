
scheduler = Scheduler();
job = submit(scheduler,{'/u/abarth/bin/run_octave.sh "svd(randn(1500,1500));"'})
isrunning(scheduler,job)
wait(scheduler,job)
disp('finished')
