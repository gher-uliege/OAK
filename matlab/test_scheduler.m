scheduler = SchedulerShell();
job = submit(scheduler,{'sleep 20'})
isrunning(scheduler,job)
wait(scheduler,job)
disp('finished')

return


scheduler = Scheduler();
job = submit(scheduler,{'/u/abarth/bin/run_octave.sh "svd(randn(1500,1500));"'})
isrunning(scheduler,job)
wait(scheduler,job)
disp('finished')
