% function stat = oak_assim_stat(initname,varnames)
% RMS and bias of assimilation increment
% Assumption: all variables in varnames have the same shape
% Input:
%   initname: init filename
%   varnames: cell of list of variables names (related to ObsXXX.names)


function stat = oak_assim_stat(initname,varnames)

init = InitFile(initname);

for v = 1:length(varnames)
  varname = varnames{v};

  rms_Hxf_yo = 0;
  bias_Hxf_yo = 0;
  rms_Hxa_yo = 0;
  bias_Hxa_yo = 0;
  stddevHxf = 0;
  stddevHxa = 0;
  
  Hxf_yo_count = 0;

  stddevHxf_yo_count = 0; 
  stddevHxa_yo_count = 0; 
 
  n = 1;

  while 1

    prefix = sprintf('Diag%03g',n);
    try
      names = get(init,sprintf('Obs%03g.names',n));
    catch      
      break
    end
    
%     if n > 5
%       warning('debuggg');       break
%     end

    i = find(strcmp(varname,names));
    
    if ~isempty(i)
      path = get(init,sprintf('Obs%03g.path',n));  
      dpath = get(init,sprintf('Diag%03g.path',n));  
      values = get(init,sprintf('Obs%03g.values',n));
      Hxf_name = get(init,sprintf('Diag%03g.Hxf',n));
      Hxa_name = get(init,sprintf('Diag%03g.Hxa',n));
      yo_name = get(init,sprintf('Diag%03g.yo',n));

      
      disp([dpath  yo_name{i}])
      
      yo  = gread([dpath  yo_name{i}]);
      Hxf = gread([dpath Hxf_name{i}]);
      Hxa = gread([dpath Hxa_name{i}]);

      if strcmp(varname,'icec')
        val = gread([path values{i}]);        
        mz = isnan(yo) & ~isnan(val);
        yo(mz) = 0;            
        Hxf(mz) = 0;            
        Hxa(mz) = 0;             
     end

     if strcmp(varname,'ui_ice') || strcmp(varname,'vi_ice')
        val = gread([path values{i}]);
        
        mz = isnan(yo) & val == 0;
        yo(mz) = NaN;
        Hxf(mz) = NaN;
        Hxa(mz) = NaN;
     end

     
     
      % forecast stat
      diff = yo - Hxf;
      Hxf_yo_count = Hxf_yo_count + ~isnan(diff);
      diff(isnan(diff)) = 0; 
      rms_Hxf_yo = rms_Hxf_yo + diff.^2;
      bias_Hxf_yo = bias_Hxf_yo + diff;

      % analysis stat
      diff = yo - Hxa;
      diff(isnan(diff)) = 0; 
      rms_Hxa_yo = rms_Hxa_yo + diff.^2;
      bias_Hxa_yo = bias_Hxa_yo + diff;       
    
      
      stddevHxf_name = get(init,sprintf('Diag%03g.stddevHxf',n));
      tmp = gread([dpath stddevHxf_name{i}]);
      if strcmp(varname,'icec')
        tmp(mz) = 0;             
      end          
      if strcmp(varname,'ui_ice') || strcmp(varname,'vi_ice')
        tmp(mz) = NaN;
      end
      
      stddevHxf_yo_count = stddevHxf_yo_count + ~isnan(tmp);
      tmp(isnan(tmp)) = 0;
      stddevHxf = stddevHxf + tmp.^2;
    
      stddevHxa_name = get(init,sprintf('Diag%03g.stddevHxa',n));
      tmp = gread([dpath stddevHxa_name{i}]);
      if strcmp(varname,'icec')
        tmp(mz) = 0;             
      end          
      if strcmp(varname,'ui_ice') || strcmp(varname,'vi_ice')
        tmp(mz) = NaN;
      end
      stddevHxa_yo_count = stddevHxa_yo_count + ~isnan(tmp);
      tmp(isnan(tmp)) = 0;
      stddevHxa = stddevHxa + tmp.^2;
    end

    n = n + 1;
  end

  rms_Hxf_yo = sqrt(rms_Hxf_yo ./ Hxf_yo_count);
  bias_Hxf_yo = bias_Hxf_yo ./ Hxf_yo_count;

  rms_Hxa_yo = sqrt(rms_Hxa_yo ./ Hxf_yo_count);
  bias_Hxa_yo = bias_Hxa_yo ./ Hxf_yo_count;

  % save in stat
  stat.(varname).('forecast').('rms') =  rms_Hxf_yo;
  stat.(varname).('forecast').('bias') =  bias_Hxf_yo;
  stat.(varname).('forecast').('stddevHx') =  sqrt(stddevHxf ./ stddevHxf_yo_count);

  stat.(varname).('analysis').('rms') =  rms_Hxa_yo;
  stat.(varname).('analysis').('bias') =  bias_Hxa_yo;
  stat.(varname).('analysis').('stddevHx') =  sqrt(stddevHxa ./ stddevHxa_yo_count);

  stat.(varname).('count') =  Hxf_yo_count;
end