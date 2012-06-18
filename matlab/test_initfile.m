initfname = '../assim.init.template';


init = InitFile(initfname);


assert(isequal(get(init,'runtype'),2))
assert(isequal(get(init,'Model.variables'),{'var1','var2'}))
assert(isequal(get(init,'Zones.corrLength.const'),[30e3,30e3]))
assert(isequal(get(init,'logfile'),'assim.log'))
%get(init,'ErrorSpace.init')

assert(isequal(get(init,'not_there',1),1))
assert(isequal(get(init,'Diag001.xf'),{'forecast.nc#var1','forecast.nc#var2'}))


assert(isequal(init.runtype,2))
assert(isequal(init.Model_variables,{'var1','var2'}))
assert(isequal(init.Zones_corrLength_const,[30e3,30e3]))
assert(isequal(init.logfile,'assim.log'))


assert(strcmp(getcopy(init,'Obs001.time'),'2010-07-06T00:30:00.00'))