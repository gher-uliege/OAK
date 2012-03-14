
currentdir = pwd;
initfile = fullfile(currentdir,'test_assim.init');

testdir = tempname;
load /home/abarth/tmp/TAvrAssim/seed.mat state
randn('state',state)

[Eic, model, fun, obs, Eforcing, t0, data, h ] = oak_create_test(testdir,initfile);


scheduler = SchedulerShell();

n = 1;
cd(testdir)
Ef = oak_assim(Eic,n,data,scheduler);
cd(currentdir)


if 0
exec = '/home/abarth/Assim/OAK/assim-gfortran-single';

n = 1;
syscmd('%s %s %d',exec,initfile,n);

path = get(init,sprintf('Diag%03g.path',n));
Eaname = get(init,sprintf('Diag%03g.Ea',n));

Ean = SVector(path,Eaname,mask,1:Nens);
end


iR = spdiag(1./(obs(1).RMSE.^2));

E = full(Eic);
for i=1:size(E,2)
  HE(:,i) = h(E(:,i));
end

[inc,info1] = ensemble_analysis(E,HE,obs(1).yo,iR);
E2 = E * info1.A;

%assert(rms(E2,full(Ef)) < 1e-5);

rms (mean(E2,2), mean(full (Ef),2))
rms(E2,full(Ef))


if 0
addpath /home/abarth/Assim/OAK/Mex/

oak_init('test_assim.init')
[E3] = oak_analysis(1,E);

rms(E2,E3)
end

