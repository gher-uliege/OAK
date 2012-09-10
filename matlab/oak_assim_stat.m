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

  n = 1;

  while 1

    prefix = sprintf('Diag%03g',n);
    try
      names = get(init,sprintf('Obs%03g.names',n));
    catch      
      break
    end
    
%     if n > 5
%       warning('debuggg')
%       break
%     end

    i = find(strcmp(varname,names));
    
    rms_Hxf_yo = 0;
    bias_Hxf_yo = 0;

    rms_Hxa_yo = 0;
    bias_Hxa_yo = 0;

    Hxf_yo_count = 0;

    if ~isempty(i)
      path = get(init,sprintf('Obs%03g.path',n));  
      dpath = get(init,sprintf('Diag%03g.path',n));  
      values = get(init,sprintf('Obs%03g.values',n));
      Hxf_name = get(init,sprintf('Diag%03g.Hxf',n));
      Hxa_name = get(init,sprintf('Diag%03g.Hxa',n));
      yo_name = get(init,sprintf('Diag%03g.yo',n));

      disp([dpath  yo_name{i}])
      
      %  yo = gread([path values{i}]);
      yo  = gread([dpath  yo_name{i}]);
      Hxf = gread([dpath Hxf_name{i}]);
      Hxa = gread([dpath Hxa_name{i}]);

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

  stat.(varname).('analysis').('rms') =  rms_Hxa_yo;
  stat.(varname).('analysis').('bias') =  bias_Hxa_yo;
end