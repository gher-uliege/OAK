
currentdir = pwd;
initfile = fullfile(currentdir,'test_assim_inflation.init');

testdir = tempname;
%testdir = [getenv('HOME') '/tmp/oak-temp'];
randn('state',0)

[t0,data,model,Eic,Eforcing, obs, fun, h ] = oak_create_test(testdir,initfile);

init = InitFile(initfile);
scale = get(init,'ErrorSpace.scale',1);

scheduler = SchedulerShell();

n = 1;
cd(testdir)
Ef = oak_assim(Eic,n,data,scheduler);
cd(currentdir)

iR = spdiag(1./(obs(1).RMSE.^2));

E = full(Eic);
for i=1:size(E,2)
  HE(:,i) = h(E(:,i));
end


E = oak_scale_ens(E,scale);
HE = oak_scale_ens(HE,scale);


[inc,info1] = ensemble_analysis(E,HE,obs(1).yo,iR);
E2 = E * info1.A;

%assert(rms(E2,full(Ef)) < 1e-5);

rms (mean(E2,2), mean(full (Ef),2))
rms (E2,full(Ef))

