function v = oak_loadens(init,key,ml)

i = find(key == '.',1,'last');
path = get(init,[key(1:i) 'path'],'');
dim = get(init,[key(1:i) 'dimension']);
names = get(init,key)
v = SVector(path,names,ml.masks,1:dim);
