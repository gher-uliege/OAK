% joins parts

function full = strjoin(parts,sep);

full = '';

for i=1:length(parts)
  if i > 1
    full = [full sep];
  end
  
  full = [full parts{i}];
end

    
