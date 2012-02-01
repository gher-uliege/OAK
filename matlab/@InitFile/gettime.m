function val = gettime(self,key,calendar)


str = get(self,key);
str = strrep(str,'T',' ');

% yyyy-mm-dd HH:MM:SS
format = 'yyyy-mm-dd HH:MM:SS';
val = datenum_cal(datevec(str,format),calendar);
