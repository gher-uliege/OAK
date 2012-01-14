function E = runoak(t0,obs,model,Eic,Eforcing)

Nens = size(Eic,2);
E = Eic;

time = [t0 time(obs)];

simulation = cell(Nens,1);

for n = 2:length(time)
  % submit ensemble run
  for i=1:Nens    
    simulation{i} = run(model,time(n-1),time(n),E(:,i),Eforcing(:,i));
  end

  % get model results
  for i=1:Nens    
    E(:,i) = result(model,simulation{i});
  end

  % assimilation
  E = oak_assim(E,n-1,obs);
end