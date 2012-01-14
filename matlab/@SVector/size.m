% extracts a column with the syntax SV(:,m)


function sz = size(self,n)

sz = [self.sv.n length(self.members)];

if nargin == 2
  sz = sz(n);
end

  