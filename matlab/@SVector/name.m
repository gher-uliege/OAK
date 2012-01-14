% name(self,i,j)
% filename of the ith variable and jth columns (ensemble member)

function s = name(self,i,j)
%keyboard
if nargin == 2
  j = 1;
end

s = fullfile(self.path,fmtprint(self.variables{i},self.members(j)));
