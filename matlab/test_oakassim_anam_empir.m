
currentdir = pwd;
initfile = fullfile(currentdir,'test_assim_anam_empir.init');

testdir = tempname;
%testdir = '/tmp/oak-temp';

randn('state',0)

[t0,data,model,Eic,Eforcing, obs, fun, h] = oak_create_test(testdir,initfile);


%% set up anamorphosis
ax = [-1e3 -1 1 1e3];
ay = [-1e4 -2 2 1e4];

init = InitFile(initfile);
transname = get(init,'Anamorphosis.transform');
path = get(init,'Anamorphosis.path');
path = fullfile(testdir,path);

gwrite(fullfile(path,transname{1}),[ax(:) ay(:)]);


fanam = @(x) interp1(ax,ay,x);
fianam = @(y) interp1(ay,ax,y);

%fanam = @(x) log(x);
%fianam = @(y) exp(y);


E = full(Eic);
Eic(:,:) = [fianam(E(1:150,:)); E(1:150,:)];

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

E(1:150,:) = fanam(E(1:150,:));
[inc,info1] = ensemble_analysis(E,HE,obs(1).yo,iR);
E2 = E * info1.A;
E2(1:150,:) = fianam(E2(1:150,:));

%assert(rms(E2,full(Ef)) < 1e-5);

rms (mean(E2,2), mean(full (Ef),2))
rms(E2,full(Ef))


if 0
addpath /home/abarth/Assim/OAK/Mex/

oak_init('test_assim.init')
[E3] = oak_analysis(1,E);

rms(E2,E3)
end

