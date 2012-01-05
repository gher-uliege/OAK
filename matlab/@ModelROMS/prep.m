% class representing a model

function prep(self,t0,t1,ic,bc,forcing,workdir)

self.workdir = workdir;

ntimes = 24*60*60*(t1 - t0)/self.dt;

[success,message,messageid] = mkdir(self.workdir);

freplace(self.template,fullfile(self.workdir,'ocean.in'), ...
        '<NtileI>',self.p.NtileI,...
        '<NtileJ>',self.p.NtileJ,...
        '<NTIMES>',ntimes,...
        '<NHIS>',self.p.nhis,...
        '<NAVG>',self.p.navg,...
        '<DT>',self.dt, ...
        '<DSTART>',t0 );


% create symbolic links for necessary files

symlink(self.p.varname,self.workdir,'delete');
symlink(self.p.grdname,self.workdir,'delete');

for i=1:length(forcing)
  symlink(forcing{i},self.workdir,'delete');
end

symlink(ic,self.workdir,'delete');
symlink(bc,self.workdir,'delete');


