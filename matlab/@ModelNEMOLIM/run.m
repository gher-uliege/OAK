% simulation = run(self,t0,t1,ic,forcing)
% make a simulation of the model self from time t0 to time t1 using
% the initial condition ic and the forcing fields forcing
%
% ic: must be a SVector (containing the initial condition). It containts EnsXXX where XXX is the 
%   ensemble member number
% forcing: must be a SVector (containing the stochastic forcing: boundary condition or atmospheric forcing)
%
% the model is run in the directory containing the initial condition
% fix me: run number is set to 48

function simulation = run(self,t0,t1,ic,forcing)


% get ensemble member number from initial condition
[S, E, TE, M, T, NM] = regexp(name(ic,1),'.*Ens([0-9]*).*');
member = str2double(T{1});

% find OPA and LIM restart file
for i=1:length(var(ic))
  str = name(ic,i);
  % tn is temperature
  if endswith(str,'tn')
    opa_restart = gfilename(str);
%  elseif endswith(str,'sxice')
%    lim_restart = gfilename(str);
  end
end
  
% working directory 
lim_restart = strrep(opa_restart,'restart','restart_ice_in');
workdir = dirname(opa_restart);

n0 = (t0 - self.torigin)*24*60*60/self.dt + 1;
n1 = (t1 - self.torigin)*24*60*60/self.dt;


opa_restart_new = strrep(opa_restart,'restart','forecast');
lim_restart_new = strrep(lim_restart,'restart','forecast');

simulation.workdir = workdir;
simulation.n0 = n0;
simulation.n1 = n1;
simulation.member = member;
olddir = pwd;

launchdir = [getenv('HOME') '/NEMO/SIMUL/ORCA2-L053-Ens'];

%cd(launchdir);


fprintf(1,'time_counter old %d\n',gread([opa_restart '#time_counter']));


args = {'srun','--nodes=1','--ntasks=1','--output',sprintf('member%03g-%d.out',member,n0)};
%keyboard
simulation.job = submit(self.scheduler,{...
                    args{:},...
                    self.script, ...
                    '--run-number',     48,...
                    '--first-step',     n0,...
                    '--last-step',      n1,...
                    '--member',         member,...
                    '--experiment',     self.experiment,...
                    '--ic-opa',         opa_restart,...
                    '--ic-lim',         lim_restart,...
                    '--out-restart-opa',opa_restart_new,...
                    '--out-restart-lim',lim_restart_new...
                   },'name',sprintf('member%03g-%d',member,n0));

%cd(olddir);
% try to prevent slurmd timeout error
system('sleep 10');

variables = var(ic);

for i = 1:length(variables)
  variables{i} = strrep(variables{i},'restart','forecast');
end

% replace the variables in SVector ic 
simulation.result = var(ic,variables);
simulation.opa_restart_new = opa_restart_new;


