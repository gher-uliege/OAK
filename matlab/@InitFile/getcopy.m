function val = getcopy(self,key,default)

ind = find(ismember(self.keys,key)==1,1,'last');

if isempty(ind)
  val = default;
else
  val = self.values{ind};
end