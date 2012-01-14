% extracts a column with the syntax SV(:,m)


function v = subsref(self,idx)

assert(strcmp(idx.type,'()'))

if strcmp(idx.subs{1},':');
  m = idx.subs{2};

  v = SVector(self.path,self.variables,self.sv.mask,self.members(m));
else
  % load all
  warning('loading all');
  tmp = full(self);
  v = subsref(tmp,idx);
end

