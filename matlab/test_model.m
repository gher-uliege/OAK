

fun = @(t0,t1,x) [x([2:end 1])];

t0 = 0;
t1 = 1;

model = ModelFun(1,fun);

x = randn(1,100);

simulation = run(model,t0,t1,x,[],[]);
xn = result(model,simulation);

xnr = fun(t0,t1,x);
assert(rms(xnr,xn) < 1e-11);

