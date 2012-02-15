fun = @(t0,t1,x,forcing) x([2:end 1],:);
%fun = @(t0,t1,x) x;

h = @(x) x([1 5],:);

model = ModelFun(1,fun);

x = randn(100,1);



Nens = 10;
time = 1:10;

E = repmat(x,[1 Nens]);
E = E + randn(size(E))/100;

obs = [];
obs(1).time = 1;
obs(1).H = h;
obs(1).yo = [1 2]';
obs(1).RMSE = [1 1]';

obs(2).time = 2;
obs(2).H = h;
obs(2).yo = [1 2]';
obs(2).RMSE = [1 1]';

t0 = 0;

Eic = E;
Eforcing = zeros(0,Nens);

data = DataSet(obs);
scheduler = [];


Ef = runoak(t0,data,model,Eic,Eforcing,scheduler);


iR = spdiag(1./(obs(1).RMSE.^2));

E = fun([],[],E);
[inc,info1] = ensemble_analysis(E,h(E),obs(1).yo,iR);
E2 = E * info1.A;

E2 = fun([],[],E2);

[inc,info2] = ensemble_analysis(E2,h(E2),obs(2).yo,iR);

E3 = E2 * info2.A;

assert(rms(E3,Ef) < 1e-5);
