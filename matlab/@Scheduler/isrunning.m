function isr = isrunning(self,job)

[status, output] = system(sprintf('qsub -j %s',job.id));

isr = status == 0;