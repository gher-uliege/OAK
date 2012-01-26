% retval = DataSetInitFile(filename,n,'properties',value,...)
% properties are:
% 'calendar': type of calendar for datenum_cal (standard, noleap,...)
% 'time_origin': time origin [year month day hour minute second]

function retval = DataSetInitFile(filename,n)

self.filename = filename;
self.init = InitFile(filename);
self.n = n;

self.exec = get(self.init,'Config.exec','assim');
self.calendar = get(self.init,'Config.calendar','standard');
str = get(self.init,'Config.time_origin','1858-11-17T00:00:00');
str = strrep(str,'T',' '); % yyyy-mm-dd HH:MM:SS  
format = 'yyyy-mm-dd HH:MM:SS';
self.t0 = datenum_cal(datevec(str,format),self.calendar);


retval = class(self,'DataSetInitFile');