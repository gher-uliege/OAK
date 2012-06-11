function val = get(self,key,default)

val=self.values{find(ismember(self.keys,key)==1,1,'last')};
