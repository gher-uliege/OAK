function v = oak_loadvec(init,key,ml)

i = find(key == '.',1,'last');
path = get(init,[key(1:i) 'path'],'');
names = get(init,key);
v = SVector(path,names,ml.masks);
