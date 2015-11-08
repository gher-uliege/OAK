
currentdir = pwd;
initfile = fullfile(currentdir,'test_assim.init');
initfile = fullfile(currentdir,'test_assim_local.init');

testdir = tempname;
%testdir = '/tmp/oak-temp';
randn('state',0)

[t0,data,model,Eic,Eforcing, obs, fun, h ] = oak_create_test(testdir,initfile);


init = InitFile(initfile);

cd(testdir)

partname = get(init,'Zones.partition');
path = get(init,'Zones.path');

modml = oak_memlayout(init,'Model');
maskv = ones(modml.VarShape(:,1));

part = gen_part(maskv);

for i=1:length(partname)
  gwrite(fullfile(testdir,path,partname{i}),part);
end


scheduler = SchedulerShell();

n = 1;
Ef = oak_assim(Eic,n,data,scheduler);

inc = full(load(init,mask(Eic),'Diag001.xa-xf'));
xf = full(load(init,mask(Eic),'Diag001.xf'));
xa = full(load(init,mask(Eic),'Diag001.xa'));
a = gread(fullfile(get(init,'Diag001.path'),get(init,'Diag001.amplitudes')));

assert( rms(inc,xa-xf) < 1e-5)
disp('increment: OK');

data.filename = fullfile(currentdir,'test_assim_local_mpi.init');

Ef2 = oak_assim(Eic,n,data,scheduler);

inc2 = full(load(init,mask(Eic),'Diag001.xa-xf'));
xf2 = full(load(init,mask(Eic),'Diag001.xf'));
xa2 = full(load(init,mask(Eic),'Diag001.xa'));
a2 = gread(fullfile(get(init,'Diag001.path'),get(init,'Diag001.amplitudes')));

cd(currentdir)


assert( rms(xf2,xf) < 1e-5)
disp('forecast: OK');

assert( rms(xa2,xa) < 1e-5)
disp('analysis: OK');

assert( rms(inc2,inc) < 1e-5)
disp('inc: OK');

assert( rms(a2,a) < 1e-5)
disp('amplitudes: OK');

assert( rms(full(Ef2),full(Ef)) < 1e-5)
disp('ensemble: OK');



