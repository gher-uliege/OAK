function val = gettime(self,key,calendar,default)

if nargin == 3
  str = get(self,key);
else
  str = get(self,key,default);
end

str = strrep(str,'T',' ');

% yyyy-mm-dd HH:MM:SS
format = 'yyyy-mm-dd HH:MM:SS';
val = datenum_cal(datevec(str,format),calendar);
