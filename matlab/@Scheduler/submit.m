function job = submit(self,args,varargin)

name = 'member';

for i=1:2:length(varargin)
  if strcmp(varargin{i},'name')
    name = varargin{i+1};
  else
    error(['unknown property: ' varargin{i}]);
  end
end

cmd = [self.command ' ' self.option_name name];

for i = 1:length(args)
  cmd = [cmd ' ' num2str(args{i})];
end

disp(cmd)
%error('stop');
[status, output] = system(cmd);

if status ~= 0
  error(['command "' cmd '" failed: ' output]);
end

if strcmp(self.command,'qsub')
  [S, E, TE, M, T]  = regexp(output,'Your job ([0-9]+) \((.*)\) has been submitted');
  job.id = T{1}{1};
  job.name = T{1}{2};
else
  [S, E, TE, M, T]  = regexp(output,'Submitted batch job ([0-9]+)');
  job.id = T{1}{1};
  job.name = name;
end

job.args = args;
job.scheduler = self;  
