%system('rm -R Ens0*');
currentdir = pwd;
initfile = fullfile(currentdir,'test_assim.init');

testdir = tempname;


init = InitFile(initfile);
[t0,data,model,Eic,Eforcing, obs, fun, h ] = oak_create_test(testdir,initfile);
E = full(Eic);

Nens = size(Eic,2);

scheduler = SchedulerShell();

cd(testdir)
Ef = runoak(t0,data,model,Eic,Eforcing,scheduler);
cd(currentdir)

iR = spdiag(1./(obs(1).RMSE.^2));

E = fun([],[],E);

for i=1:Nens
  HE(:,i) = h(E(:,i));
end

[inc,info1] = ensemble_analysis(E,HE,obs(1).yo,iR);
E2 = E * info1.A;

%rms(E2,full(Ef))
%rms (mean(E2,2), mean(full (Ef),2))
%return

E2 = fun([],[],E2);

for i=1:Nens
  HE2(:,i) = h(E2(:,i));
end

[inc,info2] = ensemble_analysis(E2,HE2,obs(2).yo,iR);

E3 = E2 * info2.A;


%assert(rms(E3,full(Ef)) < 1e-5);
rms(E3,full(Ef))
rms (mean(E3,2), mean(full (Ef),2))


