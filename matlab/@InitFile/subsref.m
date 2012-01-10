% get a value from initfile matching key
% initfile.key
% _ are subsituded by .
% e.g.
% "initfile.Model_variable" looks for Model.variable
% avoid this construct in scripts (use get(init,key) instead)

function val = subsref(self,idx)

assert(strcmp(idx.type,'.'));
key = idx.subs;
key = strrep(key,'_','.');

val = get(self,key);

