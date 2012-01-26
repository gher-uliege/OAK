
a = single(randn(100,100));
b = single(randn(100,100));

mask = ones(size(a));

path = '/tmp';
v = {'foo.nc#a','foo.nc#b'};
gwrite(fullfile(path,v{1}),a);
gwrite(fullfile(path,v{2}),b);

sv = SVector(path,v,{mask,mask});

x = full(sv);

x2 = [a(:); b(:)];
maxdiff(x2,x)


names = {'foo%03g.nc#a','foo%03g.nc#b'};
Nens = 2;

clear a b
a{1} = single(randn(100,100));
b{1} = single(randn(100,100));
a{2} = single(randn(100,100));
b{2} = single(randn(100,100));

for i=1:Nens
  gwrite(fullfile(path,sprintf(names{1},i)),a{i});
  gwrite(fullfile(path,sprintf(names{2},i)),b{i});
end

sv = SVector(path,names,{mask,mask},1:Nens);


X = full(sv);

X2(:,1) = [a{1}(:); b{1}(:)];
X2(:,2) = [a{2}(:); b{2}(:)];

maxdiff(X,X2)

maxdiff(X2(:,1),full(sv(:,1)))
maxdiff(X2(:,2),full(sv(:,2)))

c = randn(100,100);
y = [c(:); c(:)];

%sum(full(sv)(:))

sv(:,1) = y;
X2(:,1) = y;

%sum(full(sv)(:))

maxdiff(full(sv),X2)

X2(:,1:2) = [y 2*y];
sv(:,1:2) = [y 2*y];

maxdiff(full(sv),X2)


X2(:,:) = [3*y 2*y];
sv(:,:) = [3*y 2*y];

maxdiff(full(sv),X2)


A = single(randn(2,2));

maxdiff(X2 * A,sv * A)


maxdiff(X2([2:end 1],:),full(sv([2:end 1],:)))
maxdiff(X2([2:20000 1],:),full(sv([2:20000 1],:)))


names2 = {'bar%03g.nc#a','bar%03g.nc#b'};
clear a b
a{1} = single(randn(100,100));
b{1} = single(randn(100,100));
a{2} = single(randn(100,100));
b{2} = single(randn(100,100));

for i=1:Nens
  gwrite(fullfile(path,sprintf(names2{1},i)),a{i});
  gwrite(fullfile(path,sprintf(names2{2},i)),b{i});
end

sv2 = SVector(path,names2,{mask,mask},1:Nens);
sv(:,:) = sv2;

maxdiff(full(sv2),full(sv))


names3 = {'foobar%03g.nc#a','foobar%03g.nc#b'};
save(sv,path,names3);

maxdiff(gread('/tmp/foobar001.nc#a'),a{1});

%assert(isequal(var(sv),names))
%sv = var(sv,names2);
%assert(isequal(var(sv),names2))


