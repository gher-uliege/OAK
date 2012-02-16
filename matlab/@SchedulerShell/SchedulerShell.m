% job scheduler based on the Unix shell

function retval = SchedulerShell()

self = struct();
retval = class(self,'SchedulerShell',Scheduler);