% function stat = oak_assim_stat(initname,varnames)
% RMS and bias of assimilation increment
% Assumption: all variables in varnames have the same shape
% Input:
%   initname: init filename
%   varnames: cell of list of variables names (related to ObsXXX.names)


%function stat = oak_assim_stat(initname,varnames)

%end

'assim.init'
init = InitFile(initname);


for name = fieldnames(stat)'
  plotfn = plotfns.(name);
  
  plotfn(stat.(name).('forecast').('rms'),name);
  
  
end