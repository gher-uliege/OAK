function val = getcopy(self,key,default)

%problem with if key is not Obs*time
%val=self.values{find(ismember(self.keys,key)==1,1,'last')};

if nargin == 2
  val = getinitval(self.filename,key);
else
  val = getinitval(self.filename,key,default);
end

