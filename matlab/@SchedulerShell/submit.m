function job = submit(self,args,varargin)

name = 'job';

for i=1:2:length(varargin)
  if strcmp(varargin{i},'name')
    name = varargin{i+1};
  else
    error(['unknown property: ' varargin{i}]);
  end
end


cmd = '';
% expand ~ to home dir
[args{1}] = gread_tilde_expand(args{1});

for i = 1:length(args)
%    cmd = [cmd ' "' num2str(args{i}) '" ' ];
    cmd = [cmd ' ' num2str(args{i}) ' ' ];
end

cmd = [cmd ' >> ' name '.out & echo $!'];
disp(cmd);
[status,out] = system(cmd);

if status ~= 0
    error(['command "' cmd '" failed ']);
end

% str2num (instead of str2double) is necessary for octave
job.pid = str2num(out);

if isempty(job.pid)
    error(['no pid returned from  "' cmd '"']);
end
