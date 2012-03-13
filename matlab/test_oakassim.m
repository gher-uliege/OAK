
currentdir = pwd;
initfile = fullfile(currentdir,'test_assim.init');

testdir = '/tmp/tp8de9918a_1a20_4ecb_9abd_f5829d1e02b0';
%testdir = tempname;
mkdir(testdir)
cd(testdir)

init = InitFile(initfile);

sz = [10 15 1];
[x,y,z] = ndgrid(1:sz(1),1:sz(2),1:sz(3));

xname = get(init,'Model.gridX');
yname = get(init,'Model.gridY');
zname = get(init,'Model.gridZ');
path = get(init,'Model.path');
maskname = get(init,'Model.mask');

mask = {};

for i=1:length(xname)
  gwrite(fullfile(path,xname{i}),x);
  gwrite(fullfile(path,yname{i}),y);
  gwrite(fullfile(path,zname{i}),z);

  mask{i} = ones(size(x));
  gwrite(fullfile(path,maskname{i}),mask{i});
end


fun = @(t0,t1,x) x;

model = ModelFun(1,fun);


Nens = get(init,'ErrorSpace.dimension');

time = 1:2;

fmt = get(init,'ErrorSpace.init');
path = get(init,'ErrorSpace.path');

Eic = SVector(path,fmt,mask,1:Nens);

E = randn(size(Eic));
Eic(:,:) = E;


obs = [];

xobs = prctile(x(:),[40 60]);
yobs = prctile(y(:),[40 60]);
zobs = prctile(z(:),[40 60]);

h = @(v) interpn(x,y,reshape(v(1:numel(x)),size(x)),xobs,yobs);

v = gread(fullfile(path,'Ens001/model.nc#var1'));
assert(rms(h(E(:,1)),interpn(x,y,v,xobs,yobs)) < 1e-6)


for n=1:length(time)
  obs(n).time = time(n);
  obs(n).H = h;
  obs(n).yo = [1 2]';
  obs(n).RMSE = [1 1]';   

  obsprefix = sprintf('Obs%03g',n);
  obsxname = get(init,[obsprefix '.gridX']);
  obsyname = get(init,[obsprefix '.gridY']);
  obszname = get(init,[obsprefix '.gridZ']);
  obsRMSEname = get(init,[obsprefix '.rmse']);
  obsmaskname = get(init,[obsprefix '.mask']);
  obsname = get(init,[obsprefix '.value']);
  path = get(init,[obsprefix '.path']);

  gwrite(fullfile(path,obsxname{1}),xobs);
  gwrite(fullfile(path,obsyname{1}),yobs);
  gwrite(fullfile(path,obszname{1}),zobs);
  gwrite(fullfile(path,obsmaskname{1}),double(isfinite(obs(n).yo)));
  gwrite(fullfile(path,obsRMSEname{1}),obs(n).RMSE);
  gwrite(fullfile(path,obsname{1}),obs(n).yo);
end

t0 = 0;
Ebc = zeros(0,Nens);
Eforcing = zeros(0,Nens);

data = DataSetInitFile(initfile,1:length(time));

scheduler = SchedulerShell();

n = 1;
Ef = oak_assim(Eic,n,data,scheduler);



if 0
exec = '/home/abarth/Assim/OAK/assim-gfortran-single';

n = 1;
syscmd('%s %s %d',exec,initfile,n);

path = get(init,sprintf('Diag%03g.path',n));
Eaname = get(init,sprintf('Diag%03g.Ea',n));

Ean = SVector(path,Eaname,mask,1:Nens);
end





iR = spdiag(1./(obs(1).RMSE.^2));

for i=1:Nens
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

cd(currentdir)