

names = {'Ens%03g/ic.nc#x'};
names = {'ic%03g.nc#x'};
Nens = 2;
mask = ones(40,1);
path = './';

sv = SVector(path,names,{mask},1:Nens);

IC = randn(40,Nens);

sv(:,:) = IC;


rms(IC,full(sv))
