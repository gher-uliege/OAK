function xa = load(self,masks,key,members)

i = find(key == '.',1,'last');

path = realpath(get(self,[key(1:i) 'path'],'.'));
name = get(self,key);

if nargin == 3
  xa = SVector(path,name,masks);
else
  xa = SVector(path,name,masks,members);
end