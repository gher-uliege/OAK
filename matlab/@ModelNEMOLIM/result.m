% method requesting a model simulation

function result = result(self,simulation)

% wait for simulation to be finished
wait(self.scheduler,simulation.job)

% return SVector
result = simulation.result;