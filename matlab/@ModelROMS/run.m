% simulation = run(self,t0,t1,ic,forcing)
% make a model simulation
% ic: must be a SVector (containing the initial condition)
% forcing: must be a SVector (containing the stochastic forcing: boundary condition or atmospheric forcing)

% the model is run in the directory containing the initial condition

function simulation = run(self,t0,t1,ic,forcing)

str = name(ic,1);
% assume first element is zeta (surface elevation)
assert(~isempty(strfind(str,'zeta'))); 
icname = gfilename(str);
workdir = dirname(icname);

ntimes = 24*60*60*(t1 - t0)/self.dt;

[success,message,messageid] = mkdir(workdir);

freplace(self.template,fullfile(workdir,'ocean.in'), ...
        '<NtileI>',self.p.NtileI,...
        '<NtileJ>',self.p.NtileJ,...
        '<NTIMES>',ntimes,...
        '<NHIS>',self.p.nhis,...
        '<NAVG>',self.p.navg,...
        '<DT>',self.dt, ...
        '<DSTART>',t0 );


% create symbolic links for necessary files

%symlink(self.p.varname,workdir,'delete');
%symlink(self.p.grdname,workdir,'delete');

for i=1:nvars(forcing)
  frcname = gfilename(name(forcing,i));

  if ~strcmp(realpath(frcname),realpath(fullfile(workdir,basename(frcname))))
    symlink(frcname,workdir,'delete');
  end
end


if ~strcmp(realpath(icname),realpath(fullfile(workdir,'ic.nc')))
  symlink(icname,fullfile(workdir,'ic.nc'),'delete');
end


simulation.workdir = workdir;
olddir = pwd;
cd(workdir);

simulation.job = submit(self.scheduler,{self.script, 'ocean.in'});

ls('-l',workdir)
rstname = 'ocean_rst.nc';
variables = {[rstname '#zeta'],...
             [rstname '#temp'],...
             [rstname '#salt'],...
             [rstname '#u'],...
             [rstname '#v']};

% replace the variables in SVector ic 
simulation.result = var(ic,variables);
ls('-l',workdir)

cd(olddir);

