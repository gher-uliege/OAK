% retval = DataSetInitFile(filename,n,'properties',value,...)
% properties are:
% 'calendar': type of calendar for datenum_cal (standard, noleap,...)
% 'time_origin': time origin [year month day hour minute second]

function retval = DataSetInitFile(filename,n)

self.filename = filename;
self.init = InitFile(filename);
self.n = [];

if nargin == 2
  self.n = n;
end

self.exec = get(self.init,'Config.exec','assim');
self.calendar = get(self.init,'Config.calendar','standard');
self.torigin = gettime(self.init,'Config.time_origin',self.calendar,'1858-11-17T00:00:00');

retval = class(self,'DataSetInitFile');