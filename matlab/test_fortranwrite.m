
fmt = '("Ens",I3.3,"/model.nc#var1")';

assert(isequal(fortranwrite(fmt,1),'Ens001/model.nc#var1'));
assert(isequal(fmtprint(fmt,1),'Ens001/model.nc#var1'));
assert(isequal(fmtprint('Ens%03g/model.nc#var1',1),'Ens001/model.nc#var1'));