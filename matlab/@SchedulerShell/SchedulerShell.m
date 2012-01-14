% class integrates with the SGE scheduler

function retval = SchedulerShell()

self = struct();
retval = class(self,'SchedulerShell',Scheduler);