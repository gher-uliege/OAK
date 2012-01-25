
cd /u/abarth/tmp/NEMO/ORCA2/Ensemble-test_1p8_20_70_short

experiment = 'Ensemble-test_1p8_20_70_short';


%'gher-diva:/home/abarth/tmp/NEMO/ORCA2/OUTPUT/ORCA2-L053/restart/restart.ORCA2-L053.P0047.nc'
%'gher-diva:/home/abarth/tmp/NEMO/ORCA2/OUTPUT/ORCA2-L053/restart/restart_ice_in.ORCA2-L053.P0047.nc'

cal='noleap';
dt = 5760;
torigin = datenum_cal(1960,01,01,0,0,0,cal);
t0 = datenum_cal(2007,1,1,0,0,0,cal);
t1 = datenum_cal(2008,1,1,0,0,0,cal);

t1 = datenum_cal(2007,1,3,0,0,0,cal);


scheduler = Scheduler();
script = '/u/abarth/NEMO/SIMUL/ORCA2-L053-Ens/run-ens3.bash ';
%script = 'echo ';


Nens = 3;
initfile = '/u/abarth/tmp/NEMO/ORCA2/Ensemble-test_1p8_20_70_short/assim.init';
data = DataSetInitFile(initfile,1);
init = InitFile(initfile);


mpath = get(init,'Model.path');
maskname = get(init,'Model.mask');
masks = {};
for i=1:length(maskname)
  masks{i} = gread(fullfile(mpath,maskname{i}));
end

  
path = get(init,'Diag001.path');
path = realpath(path);
fmt = get(init,'Diag001.Ea');


Eic = SVector(path,fmt,masks,1:Nens);

model = ModelNEMOLIM(dt,script,experiment,torigin,cal,scheduler);

Eforcing = zeros(0,Nens);

simulation = {};

for i=1:Nens
    simulation{i} = run(model,t0,t1,Eic(:,i),Eforcing(:,i));
end

%xn = result(model,simulation);