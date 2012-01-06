function isr = isrunning(self,job)

%[status, output] = system(sprintf('qsub -j %s',job.id));

[status, output] = system(sprintf('qstat | awk ''$1 == %s { print $5 }''',job.id)); 

isr = ~isempty(findstr(output,'r')) || ~isempty(findstr(output,'qw'));