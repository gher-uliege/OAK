function E = runoak(t0,obs,model,Eic,Eforcing,scheduler)

Nens = size(Eic,2);
E = Eic;

timeobs = [t0 time(obs)];

simulation = cell(Nens,1);

% time loop
for n = 2:length(timeobs)
  % submit ensemble run
  for i=1:Nens    
    simulation{i} = run(model,timeobs(n-1),timeobs(n),E(:,i),Eforcing(:,i));
  end

  % get model results
  for i=1:Nens    
    E(:,i) = result(model,simulation{i});
  end

  % assimilation
  E = oak_assim(E,n-1,obs,scheduler);
end