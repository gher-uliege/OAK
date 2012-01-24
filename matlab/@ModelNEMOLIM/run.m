% simulation = run(self,t0,t1,ic,forcing)
% make a model simulation
% ic: must be a SVector (containing the initial condition)
% forcing: must be a SVector (containing the stochastic forcing: boundary condition or atmospheric forcing)

% the model is run in the directory containing the initial condition

function simulation = run(self,t0,t1,ic,forcing)

% find OPA and LIM restart file
for i=1:length(var(ic))
  str = name(ic,i);
  % tn is temperature
  if endswith(str,'tn')
    opa_restart = gfilename(str);
  elseif endswith(str,'sxice')
    lim_restart = gfilename(str);
  end
end
  
% fix me ? which workdir?
workdir = dirname(opa_restart);

ntimes = 24*60*60*(t1 - t0)/self.dt;


n0 = (t0 - self.torigin)*24*60*60/self.dt + 1
n1 = (t1 - self.torigin)*24*60*60/self.dt 



% create symbolic links for necessary files

%symlink(self.p.varname,workdir,'delete');
%symlink(self.p.grdname,workdir,'delete');



simulation.workdir = workdir;
olddir = pwd;
cd(workdir);

simulation.job = submit(self.scheduler,{self.script, 
                    '--run-number',48,...
                    '--first-step',n0,...
                    '--last-step', n1,...
                    '--member', 102,...
                    '--experiment', self.experiment,...
                    '--ic-opa', opa_restart,...
                    '--ic-lim', lim_restart,...
                    '--out-restart-opa','~/tmp/NEMO/ORCA2/Ensemble-test_1p8_20_70_short/Ens102/restart.nc',...
                    '--out-restart-lim','~/tmp/NEMO/ORCA2/Ensemble-test_1p8_20_70_short/Ens102/restart_ice_in.nc'...
                   });    

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

