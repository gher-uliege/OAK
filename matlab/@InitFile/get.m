% val = get(init,key)
% val = get(init,key,default)
%
% get and parses the value corresponding to key in the initfile.
% An error is produced is key does not exist unless the default value
% is specified.

function val = get(self,key,default)

%{
if nargin == 2
  val = getinitval(self.filename,key);
else
  val = getinitval(self.filename,key,default);
end

%}

% loop over all keys in initfile

found = false;

for i=1:length(self.keys2) 
    if fnmatch(self.keys2{i},key)
        val = parseVal(self,self.values2{i});
        found = true;
    end
end

if ~found
  if nargin == 3
    val = default;
  else    
    error(['key ' key ' not found in ' initfile]);
  end
end
        
