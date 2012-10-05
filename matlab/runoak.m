% E = runoak(t0,obs,model,Eic,Eforcing,scheduler)
%
% Run an ensemble simulation and assimilate the observation when they are available (determined by
% time(obs))
%
% Input:
%
% t0: start time (calendar and time origin are given in the init file)
% obs: observations (instance of DataSetInitFile or DataSet)
% model: model object implementing methods run and result
% Eic: ensemble of initial conditions of the size n x Nens where n is the state vector size
%   and Nens the number of ensemble members. Eic is either a matlab array or an instance of SVector
% Eforcing: ensemble of forcing fields of the size p x Nens where p is the size of a forcing vector.
%   Eforcing is either a matlab array or an instance of SVector. p can be zero no forcings are perturbed.
% scheduler: instance of Scheduler or SchedulerShell used to run ensemble simulation and assimilation on
%   multiple processors.
%
% Output:
% E: ensemble of states at final time
%
% Notes:
% more information is available at 
% http://modb.oce.ulg.ac.be/mediawiki/index.php/Ocean_Assimilation_Kit

% Author Alexander Barth <a.barth@ulg.ac.be> 2012

function E = runoak(t0,obs,model,Eic,Eforcing,scheduler)

Nens = size(Eic,2);
E = Eic;

timeobs = [t0 time(obs)];

simulation = cell(Nens,1);

% time loop
for n = 2:length(timeobs)
  
  if timeobs(n-1) ~= timeobs(n)
    % submit ensemble run
    for i=1:Nens    
      simulation{i} = run(model,timeobs(n-1),timeobs(n),E(:,i),Eforcing(:,i));
    end

    % get model results
    for i=1:Nens    
      E(:,i) = result(model,simulation{i});
    end
  end
  
  % assimilation
  E = oak_assim(E,n-1,obs,scheduler);
end

% LocalWords:  runoak Eic Eforcing init DataSetInitFile DataSet Nens matlab SVector SchedulerShell
