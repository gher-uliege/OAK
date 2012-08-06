% E = oak_scale_ens(E,scale)
% inflate ensemble E by factor scale

function E = oak_scale_ens(E,scale)

meanE = mean(E,2);

for i=1:size(E,2)
    E(:,i) = meanE + scale * (E(:,i) - meanE);
end
