% la = oak_memlayout(init,key)
% la = oak_memlayout({mask1,mask2,...})
% load  memory layout based on key from InitFile init

function la = oak_memlayout(init,key)

if nargin == 2
  path = get(init,[key '.path']);
  masksnames = get(init,[key '.mask']);

  masks = cell(1,length(masksnames));
  for m = 1:length(masksnames)
    masks{m} = ~isnan(gread(fullfile(path,masksnames{m})));
  end
else
  masks = init;
end


la.StartIndex(1) = 1;
la.StartIndexSea(1) = 1;
la.nvar = length(masks);
la.VarShape = ones(5,la.nvar);


for m = 1:length(masks)
  %la.VarShape{m} = size(masks{m});
  sz = size(masks{m});
  la.VarShape(1:length(sz),m) = sz;
  la.ndim = 1:length(sz);
  
  la.VarSize(m) = numel(masks{m});
  la.VarSizeSea(m) = sum(masks{m}(:) == 1);    
  
  la.EndIndex(m) = la.StartIndex(m) + la.VarSize(m)-1;
  la.EndIndexSea(m) = la.StartIndexSea(m) + la.VarSizeSea(m)-1;

  if m ~= length(masks) 
    la.StartIndex(m+1) = la.EndIndex(m)+1;
    la.StartIndexSea(m+1) = la.EndIndexSea(m)+1;
  end
end

la.TotSizeSea = sum(la.VarSizeSea);
la.TotSize = sum(la.VarSize);
la.effsize = la.TotSizeSea;

for m=1:length(masks)
  la.mask(la.StartIndex(m):la.EndIndex(m)) = masks{m}(:);
end

la.SeaIndex = -ones(la.TotSize,1);

ind = find(la.mask);
la.SeaIndex(ind) = 1:length(ind);
la.InvIndex = ind;

%j = 1;
%for i=1:la.TotSize
%  if la.mask(i) == 1
%    la.SeaIndex(i) = j;
%    la.InvIndex(j) = i;
%    j=j+1;
%  end
%end

la.masks = masks;