function retval = ModelFun(dt,fun)

self.dt = dt;
self.fun = fun;

retval = class(self,'ModelFun',Model(dt));