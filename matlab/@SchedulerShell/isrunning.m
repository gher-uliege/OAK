function isr = isrunning(self,job)

cmd = sprintf('ps %d',job.pid);
[ret,out] = system(cmd);

isr = ret == 0;

