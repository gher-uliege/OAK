% class representing a ROMS model

% template
% workdir

% param should contain:
% NtileI
% NtileJ
% ntimes
% nhis
% navg
% varname (string)
% grdname (string)
% frcname (cell array)

function retval = ModelROMS(dt,template,param,scheduler)

self.dt = dt;
self.template = template;
self.p = param;
self.scheduler = scheduler;

retval = class(self,'ModelROMS',Model(dt));