% class representing a model

function simulation = prep(self,t0,t1,ic,bc,forcing,workdir)

simulation.model = self;
simulation.workdir = workdir;

ntimes = 24*60*60*(t1 - t0)/self.dt;

[success,message,messageid] = mkdir(workdir);

freplace(self.template,fullfile(workdir,'ocean.in'), ...
        '<NtileI>',self.p.NtileI,...
        '<NtileJ>',self.p.NtileJ,...
        '<NTIMES>',ntimes,...
        '<NHIS>',self.p.nhis,...
        '<NAVG>',self.p.navg,...
        '<DT>',self.dt, ...
        '<DSTART>',t0 );


% create symbolic links for necessary files

symlink(self.p.varname,workdir,'delete');
symlink(self.p.grdname,workdir,'delete');

for i=1:length(forcing)
  symlink(forcing{i},workdir,'delete');
end

symlink(ic,workdir,'delete');
symlink(bc,workdir,'delete');


