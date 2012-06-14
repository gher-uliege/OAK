% class integrates with the SGE scheduler

function retval = Scheduler(type)

if nargin == 0
  type = 'SGE';
end

self.type = type;

if strcmp(type,'SGE')
  self.command = 'qsub';
  self.option_name = '-N ';
elseif strcmp(type,'SLURM')
  self.command = 'sbatch';
  self.option_name = '--job-name=';
else
  warning(['unknown type ' type]);
end


retval = class(self,'Scheduler');