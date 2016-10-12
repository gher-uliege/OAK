
currentdir = pwd;
initfile = fullfile(currentdir,'test_assim.init');
%initfile = fullfile(currentdir,'test_assim_local.init');

testdir = tempname;
%testdir = [getenv('HOME') '/tmp/oak-temp'];
randn('state',0)

[t0,data,model,Eic,Eforcing, obs, fun, h] = oak_create_test(testdir,initfile);


scheduler = SchedulerShell();

n = 1;
cd(testdir)
% run OAK
Ef = oak_assim(Eic,n,data,scheduler);

% load diagnostics
init = InitFile(initfile);
obsml = oak_memlayout(init,'Obs001');
modml = oak_memlayout(init,'Model');

stddevxf = full(oak_loadvec(init,['Diag001.stddevxf'],modml));
stddevxa = full(oak_loadvec(init,['Diag001.stddevxa'],modml));

stddevHxf = full(oak_loadvec(init,['Diag001.stddevHxf'],obsml));
stddevHxa = full(oak_loadvec(init,['Diag001.stddevHxa'],obsml));
cd(currentdir)

iR = spdiag(1./(obs(1).RMSE.^2));

E = full(Eic);
for i=1:size(E,2)
  HE(:,i) = h(E(:,i));
  HEa(:,i) = h(Ef(:,i));
end

[inc,info1] = ensemble_analysis(E,HE,obs(1).yo,iR);
E2 = E * info1.A;

%assert(rms(E2,full(Ef)) < 1e-5);

rms (mean(E2,2), mean(full (Ef),2))
rms(E2,full(Ef))

for i=1:size(E2,2)
  HEa(:,i) = h(E2(:,i));
end

% check diagnostics


assert( rms(stddevxf,std(E,[],2)) < 1e-5)
disp('stddevxf: OK');

assert( rms(stddevxa,std(E2,[],2)) < 1e-5)
disp('stddevxa: OK');

assert( rms(stddevHxf,std(HE,[],2)) < 1e-5)
disp('stddevHxf: OK');

assert( rms(stddevHxa,std(HEa,[],2)) < 1e-5)
disp('stddevHxa: OK');