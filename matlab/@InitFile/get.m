function val = get(self,key,default)

if nargin == 2
  val = getinitval(self.filename,key);
else
  val = getinitval(self.filename,key,default);
end