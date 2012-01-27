% method requesting a model simulation

function result = result(self,simulation)

% wait for simulation to be finished
wait(self.scheduler,simulation.job)

% return SVector
result = simulation.result;

% check time counter
time_counter = gread([simulation.opa_restart_new '#time_counter']);

if time_counter ~= simulation.n1
    error('time counter is %d while it should be %d\n',time_counter,simulation.n1);
end
