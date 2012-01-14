function v = var(self,newv)

if nargin == 1
  v = self.variables;
else
  self.variables = newv;
  v = self;
end
