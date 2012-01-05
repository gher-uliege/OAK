% class integrates with the SGE sheduler

function retval = Job()


self.command = 'qsub';

retval = class(self,'Job');