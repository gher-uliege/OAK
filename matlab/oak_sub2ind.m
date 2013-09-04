function [index,valid] = oak_sub2ind(ML,v,i,j,k,valid) 

index = -1;
valid = 1 <= v & v <= ML.nvar;

if valid
  valid = 1 <= i & i <= ML.VarShape(1,v) &  ...
          1 <= j & j <= ML.VarShape(2,v) &  ...
          1 <= k & k <= ML.VarShape(3,v);
  
  if valid
    index = ML.StartIndex(v) + i-1 + ML.VarShape(1,v) * (j-1 + ML.VarShape(2,v) * (k-1));
    index = ML.SeaIndex(index);
    valid =  index ~= -1;
  end
end
  
