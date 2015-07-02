
n = 8;
Nens = 100;
fv = -9999999;
fname = 'domain.nc';

delete(fname)
nccreate(fname,'mask','Dimensions',{'x',n},'FillValue',fv);
nccreate(fname,'lon','Dimensions',{'x'},'FillValue',fv);
nccreate(fname,'lat','Dimensions',{'x'},'FillValue',fv);
nccreate(fname,'z','Dimensions',{'x'},'FillValue',fv);
nccreate(fname,'part','Dimensions',{'x'},'FillValue',fv);
nccreate(fname,'ens','Dimensions',{'x',n,'ensemble',Nens},'FillValue',fv);


ncwrite(fname,'mask',ones(n,1))
ncwrite(fname,'lon',[1:n]')
ncwrite(fname,'lat',zeros(n,1))
ncwrite(fname,'z',zeros(n,1))
ncwrite(fname,'part',[1:n]')

E = repmat([1:n]',[1 Nens]) + 0.01 * randn(n,Nens);
ncwrite(fname,'ens',E)


obsfname = 'obs1.nc';
m = 1;
delete(obsfname)
nccreate(obsfname,'mask','Dimensions',{'x',m},'FillValue',fv);
nccreate(obsfname,'lon','Dimensions',{'x'},'FillValue',fv);
nccreate(obsfname,'lat','Dimensions',{'x'},'FillValue',fv);
nccreate(obsfname,'z','Dimensions',{'x'},'FillValue',fv);
nccreate(obsfname,'rmse','Dimensions',{'x'},'FillValue',fv);
nccreate(obsfname,'var1','Dimensions',{'x'},'FillValue',fv);


ncwrite(obsfname,'mask',ones(m,1))
ncwrite(obsfname,'lon',[4]')
ncwrite(obsfname,'lat',zeros(m,1))
ncwrite(obsfname,'z',zeros(m,1))
ncwrite(obsfname,'rmse',1)
ncwrite(obsfname,'var1',1)
