function E = oak_assim(E,n,obs,scheduler);

fprintf(1,'Assimilation %d \n',n);
Nens = size(E,2);

if ~isa(obs,'DataSetInitFile')
  
  yo = obs(n).yo;
  iR = spdiag(1./(obs(n).RMSE.^2));
  
  % extract observed values from ensemble
  
  for j=1:Nens
    Hx = obs(n).H(E(:,j));
    
    if (j==1)
      HE = zeros([numel(Hx) Nens]);
    end
    
    HE(:,j) = Hx;
  end
  
  %keyboard
  [inc,info] = ensemble_analysis(full(E),HE,yo,iR);
  E(:,:) = E * info.A;
else    
  initfile = obs.filename;       
  init = InitFile(initfile);
  
  
  if get(init,'Config.coupling',1) == 1
    % coupling with shell
    
    masks = mask(E);
    
    % forecast ensemble
    Efname = get(init,'ErrorSpace.init');
    Efpath = get(init,'ErrorSpace.path');
    Ef = SVector(Efpath,Efname,masks,1:Nens);
    
    % make current ensemble available at Efpath,Efname
    save(E,Efpath,Efname);
    %Ef(:,:) = E;
    
    % analysis ensemble
    Eapath = realpath(get(init,sprintf('Diag%03g.path',n)));
    Eaname = get(init,sprintf('Diag%03g.Ea',n));
    Ea = SVector(Eapath,Eaname,masks,1:Nens);
    
    % copy current ensemble Efname over to Eaname
    % so that all other variables which should not be modified by the
    % assimilation are there
    save(Ef,Eapath,Eaname,'copy');
    exec = get(init,'Config.exec');
    
    job = submit(scheduler,{exec,...
                        initfile,n},'name',sprintf('analysis%03g',n));
    wait(scheduler,job)
    
    %syscmd('%s %s %d',exec,initfile,n);
    
    E = Ea;
  else
    % coupling with mex

    oak_init(initfile);
    [E] = oak_analysis(n,full(E));
    oak_done();
  end
  
end