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

function retval = ModelROMS(dt,template,param)

self.dt = dt;
self.template = template;
self.p = param;
self.workdir = '';

retval = class(self,'ModelROMS',Model(dt));