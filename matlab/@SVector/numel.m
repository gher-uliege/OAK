% number of elements 

function n = numel(self)

n = self.sv.n * length(self.members);
