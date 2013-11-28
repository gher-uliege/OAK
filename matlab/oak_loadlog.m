function res = oak_loadlog(filename,var,step,stat);
%function res = oak_loadlog(filename,pat);

% fname = [tempname '.txt'];
% %pat = '^...\.forecast.rms_yo-Hx';
% cmd = sprintf('awk ''/%s/ { print $2 } '' "%s" > %s',pat,filename,fname);

% disp(cmd);
% [o,ret] = system(cmd);

% res = load(fname);

% delete(fname);


% filename = 'assim.log-00001';

% tline = '001.temp.forecast.bias_yo-Hx                -0.9527844E-02'
% [S, E, TE, M, T, NM] = regexp(tline,'(?<index>\d{3,3})\.(?<var>\w*)\.(?<step>\w*)\.(?<stat>[^. ]*) +(?<val>.*)');

% var = 'temp';
% step = 'forecast';
% stat = 'rms_yo-Hx';

mask = [];
fid=fopen(filename,'r');
while 1
    tline = fgetl(fid);
    if ~ischar(tline),
      
       break,
    else
       [S, E, TE, M, T, NM] = regexp(tline,'(?<index>\d+)\.(?<var>\w*)\.(?<step>\w*)\.(?<stat>[^. ]*) +(?<val>.*)');

       if S
         if strcmp(var,NM.var) && strcmp(step,NM.step) && strcmp(stat,NM.stat) 
           res(str2double(NM.index)) = str2double(NM.val);
           mask(str2double(NM.index)) = 1;
         end
       end
    end
end
fclose(fid);
res(mask == 0) = NaN;