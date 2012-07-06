% ml = oak_memlayout(init,key)
% load  memory layout based on key from InitFile init

function ml = oak_memlayout(init,key)

path = get(init,[key '.path'])
masksnames = get(init,[key '.mask']);

masks = cell(1,length(masksnames));

for i=1:length(masksnames)
  masks{i} = ~isnan(gread(fullfile(path,masksnames{i})));
end

ml.masks = masks;