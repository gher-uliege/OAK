function [t0,data,model,Eic,Eforcing, obs, fun, h] = oak_create_test(testdir,initfile,sz)


if nargin == 2
%  sz = [10 15 3];
  sz = [10 15 1];
end

ndim = ndims(ones(sz));

[success,message] = mkdir(testdir);
[success,message] = mkdir(fullfile(testdir,'Common'));
[success,message] = mkdir(fullfile(testdir,'Obs'));
[success,message] = mkdir(fullfile(testdir,'Analysis001'));
[success,message] = mkdir(fullfile(testdir,'Analysis002'));


init = InitFile(initfile);

[x,y,z] = ndgrid(1:sz(1),1:sz(2),1:sz(3));

xname = get(init,'Model.gridX');
yname = get(init,'Model.gridY');
zname = get(init,'Model.gridZ');
path = get(init,'Model.path');
maskname = get(init,'Model.mask');

path = fullfile(testdir,path);
mask = {};

for i=1:length(xname)
  gwrite(fullfile(path,xname{i}),x);
  gwrite(fullfile(path,yname{i}),y);
  gwrite(fullfile(path,zname{i}),z);

  mask{i} = ones(size(x));
  gwrite(fullfile(path,maskname{i}),mask{i});
end


fun = @(t0,t1,x,forcing) x;

model = ModelFun(1,fun);


Nens = get(init,'ErrorSpace.dimension');

time = 1:2;

fmt = get(init,'ErrorSpace.init');
path = get(init,'ErrorSpace.path');
path = fullfile(testdir,path);

Eic = SVector(path,fmt,mask,1:Nens);
%full(Eic)(1,1)

E = randn(size(Eic));
Eic(:,:) = E;
assert(rms(E,full(Eic)) < 1e-5)

obs = [];

xobs = prctile(x(:),[40 60]);
yobs = prctile(y(:),[40 60]);
zobs = prctile(z(:),[40 60]);


v = gread(fullfile(testdir,'Ens001/model.nc#var1'));

if ndim == 2
  h = @(v) interpn(x,y,reshape(v(1:numel(x)),size(x)),xobs,yobs);
  assert(rms(h(E(:,1)),interpn(x,y,v,xobs,yobs)) < 1e-6)
else
  h = @(v) interpn(x,y,z,reshape(v(1:numel(x)),size(x)),xobs,yobs,zobs);
  assert(rms(h(E(:,1)),interpn(x,y,z,v,xobs,yobs,zobs)) < 1e-6)
end


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
  path = fullfile(testdir,path);

  gwrite(fullfile(path,obsxname{1}),xobs);
  gwrite(fullfile(path,obsyname{1}),yobs);
  gwrite(fullfile(path,obszname{1}),zobs);
  gwrite(fullfile(path,obsmaskname{1}),double(isfinite(obs(n).yo)));
  gwrite(fullfile(path,obsRMSEname{1}),obs(n).RMSE);
  gwrite(fullfile(path,obsname{1}),obs(n).yo);
end


t0 = 0;
Eforcing = zeros(0,Nens);

data = DataSetInitFile(initfile,1:length(time));
