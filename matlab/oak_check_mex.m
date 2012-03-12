addpath /home/abarth/Assim/OAK/Mex/

%oakmex('init')
%return

oak_init('../matlab/test_assim.init')

n = 10*15*2;
r = 3;

Ef = reshape(1:n*r,n,r);

[Ea] = oak_analysis(1,Ef);

Ea

oak_done