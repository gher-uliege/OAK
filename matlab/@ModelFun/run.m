% method requesting a model simulation

function simulation = run(self,t0,t1,ic,forcing)

simulation.result = self.fun(t0,t1,ic,forcing);
