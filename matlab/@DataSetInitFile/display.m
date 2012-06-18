% retval = DataSetInitFile(filename,n,'properties',value,...)
% properties are:
% 'calendar': type of calendar for datenum_cal (standard, noleap,...)
% 'time_origin': time origin [year month day hour minute second]

function display(self)

fprintf(1,'Calendar: %s\n',self.calendar);
fprintf(1,'torigin:  %d\n',self.torigin);

