function E = oak_ens(t0,t1,model,Eic,Eforcing,scheduler)

% submit ensemble run
for i=1:Nens    
  simulation{i} = run(model,t0,t1,E(:,i),Eforcing(:,i));
end

% get model results
for i=1:Nens    
  E(:,i) = result(model,simulation{i});
end
