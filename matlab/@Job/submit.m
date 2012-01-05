function j = submit(self,args)


cmd = self.command;

for i = 1:length(args)
  cmd = [cmd ' ' args{i}];
end

[status, output] = system(cmd);

if status ~= 0
  error(['command "' cmd '" failed: ' output]);
end

[S, E, TE, M, T]  = regexp(output,'Your job ([0-9]+) \((.*)\) has been submitted');

job.id = T{1}{1};
job.name = T{1}{2};
job.args = args;