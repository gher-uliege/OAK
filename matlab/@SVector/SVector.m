function retval = SVector(path,variables,masks,members)

self.path = path;
self.variables = variables;
self.sv = [];
self.members = [];
self.masks = [];

if nargin >= 3
  self.masks = masks;
  self.sv = statevector_init(masks{:});
end

if nargin >= 4
  self.members = members;
end

retval = class(self,'SVector');