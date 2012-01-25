% retval = SVector(path,variables)
% retval = SVector(path,variables,masks)
% retval = SVector(path,variables,masks,members)
%
% creates a SVector object
%
% Input:
%   masks: cell array of masks (or filenames of masks)

function retval = SVector(path,variables,masks,members)

self.path = path;
self.variables = variables;
self.sv = [];
self.members = [];
self.masks = [];

if nargin >= 3
  
  if ischar(masks{1})
    % load mask first if filenames are given
    names = masks;
    for i=1:length(names)
      masks{i} = gread(names{i});
    end
  end
  
  self.masks = masks;
  self.sv = statevector_init(masks{:});
end

if nargin >= 4
  self.members = members;
end

retval = class(self,'SVector');