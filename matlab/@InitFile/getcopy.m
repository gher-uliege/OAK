function val = getcopy(self,key,default)

ind = find(ismember(self.keys,key)==1,1,'last');

if length(ind) == 0
  val = default;
else
  val = self.values{ind};
end