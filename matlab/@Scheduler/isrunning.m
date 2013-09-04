function isr = isrunning(self,job)

%[status, output] = system(sprintf('qsub -j %s',job.id));

if strcmp(self.command,'qsub')
  [status, output] = system(sprintf('qstat | awk ''$1 == %s { print $5 }''',job.id)); 
  isr = ~isempty(findstr(output,'r')) || ~isempty(findstr(output,'qw')) || ~isempty(findstr(output,'t'));  
else
  [status, output] = system(sprintf('squeue --jobs=%s --format=%%T',job.id)); 
  isr = ~isempty(findstr(output,'RUNNING')) ||  ~isempty(findstr(output,'PENDING'));  
end
