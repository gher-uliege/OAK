% class representing a model

function job = run(self,scheduler)

olddir = pwd;
cd(self.workdir);

job = submit(scheduler,{self.script, 'ocean.in'});

cd(olddir);

