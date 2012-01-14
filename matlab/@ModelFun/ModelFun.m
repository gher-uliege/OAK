% Create ModelFun object for a model which is available as a callback
% function fun
% dt: time step
% result = self.fun(t0,t1,ic,forcing);

function retval = ModelFun(dt,fun)

self.dt = dt;
self.fun = fun;

retval = class(self,'ModelFun',Model(dt));