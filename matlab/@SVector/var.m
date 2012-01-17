function v = var(self,newv)

if nargin == 1
  v = self.variables;
else
  v = self;
  v.variables = newv;
end
