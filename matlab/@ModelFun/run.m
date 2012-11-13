% method requesting a model simulation

function simulation = run(self,t0,t1,ic,forcing)

fprintf('Model integration from %g to %g \n',t0,t1);
simulation.result = self.fun(t0,t1,ic,forcing);
