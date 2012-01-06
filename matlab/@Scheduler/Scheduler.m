% class integrates with the SGE scheduler

function retval = Scheduler()


self.command = 'qsub';

retval = class(self,'Scheduler');