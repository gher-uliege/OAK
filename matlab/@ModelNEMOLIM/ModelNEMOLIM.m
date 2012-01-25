% class representing a ROMS model
%
% template
% workdir
%
% param should contain:
% NtileI
% NtileJ
% nhis
% navg

%% not any more
%% Varname (string)
%% grdname (string)
%% frcname (cell array)

function retval = ModelNEMOLIM(dt,script,experiment,torigin,cal,scheduler)

self.dt = dt;
self.script = script;
self.experiment = experiment;
self.torigin = torigin;
self.cal = cal;
self.scheduler = scheduler;

retval = class(self,'ModelNEMOLIM',Model(dt));