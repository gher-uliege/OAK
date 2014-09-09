
currentdir = pwd;
initfile = fullfile(currentdir,'test_assim.init');
initfile = fullfile(currentdir,'test_assim_local.init');

testdir = tempname;
testdir = '/tmp/oak-temp';
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

a = gread(fullfile(get(init,'Diag001.path'),get(init,'Diag001.amplitudes')));

data.filename = fullfile(currentdir,'test_assim_local_mpi.init');

Ef2 = oak_assim(Eic,n,data,scheduler);
a2 = gread(fullfile(get(init,'Diag001.path'),get(init,'Diag001.amplitudes')));

cd(currentdir)


assert( rms(a2,a) < 1e-5)
disp('amplitudes: OK');

assert( rms(full(Ef2),full(Ef)) < 1e-5)
disp('ensemble: OK');
