function job = submit(self,args)


cmd = '';

for i = 1:length(args)
  cmd = [cmd ' ' args{i}];
end


cmd = [cmd ' > /dev/null & echo $!'];
%disp(cmd);
[status,out] = system(cmd);

job.pid = str2num(out);

if status ~= 0
  error(['command "' cmd '" failed ']);
end
