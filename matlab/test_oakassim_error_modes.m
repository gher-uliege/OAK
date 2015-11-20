clear mask

currentdir = pwd;
initfile = fullfile(currentdir,'test_assim_error_modes.init');

testdir = tempname;
%testdir = '/tmp/oak-temp';

randn('state',0)

[t0,data,model,Eic,Eforcing, obs, fun, h] = oak_create_test(testdir,initfile);

init = InitFile(initfile);

scale = get(init,'ErrorSpace.scale',1);

scheduler = SchedulerShell();

n = 1;
cd(testdir)
delete('analysis001.out')
Ef = oak_assim(Eic,n,data,scheduler);

xa2 = load(init,mask(Eic),sprintf('Diag%03g.xa',n));
Sa2 = load(init,mask(Eic),sprintf('Diag%03g.Sa',n),1:size(Eic,2));
stddevxf = load(init,mask(Eic),sprintf('Diag%03g.stddevxf',n));

type('analysis001.out')

cd(currentdir)

Sa2 = full(Sa2);
Pa2 = Sa2*Sa2';

R = spdiag(obs(1).RMSE.^2);

n = size(Eic,1);
xf = 2*ones(n,1);
S = scale * full(Eic);
Pf = S*S';

I = eye(n,n);
yo = obs(1).yo;
H = zeros(length(yo),n);

for i=1:n
  H(:,i) = h(I(:,i));
end

[xa,Pa] = simple_assim(xf,Pf,yo,R,H);

rms (full(xa2), xa)
rms(Pa2,Pa)
rms(full(stddevxf),sqrt(diag(Pf)))
