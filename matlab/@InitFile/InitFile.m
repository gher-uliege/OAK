% init = InitFile(filename)
% Class representing a initilization file, i.e. a list of key-value pairs
%


function retval = InitFile(filename)

key='Obs.\w*.time';

self.filename = filename;

% all keys and values
self.all_keys = {};
self.all_values = {};

% keys and values matching key
self.keys = {};
self.values = {};


% open the file

fid=fopen(filename,'r');
if fid == -1
    %error opening file
    error(['unable to open ' filename]);
end

tline=fgetl(fid);

while ischar(tline)  %Read all the file
    tline=fgetl(fid);
    if ~ischar(tline)
        break;
    end
    
    line = strtrim(tline);
    
    if isempty(line)
        continue
    elseif line(1) == '!' || line(1) == '#'
        continue
    else
        p = strsplit(line,'=');
        k = strtrim(p{1});
        v = strtrim(p{2});
        
        if length(p) == 2
            self.all_keys{end+1} = k;
            self.all_values{end+1} = v;
        end
        
        if ~isempty(regexp(tline,key,'match'))
            self.keys{end+1} = k;
            self.values{end+1} = parseVal_(v);
        end
    end
end

fclose(fid);

retval = class(self,'InitFile');

end

function val = parseVal_(str)

val = [];
if str(1) == '[' && str(end) == ']'
    str(1) = '{';
    str(end) = '}';
end

eval(['val = ' str ';']);
if iscell(val)
    if isnumeric(val{1})
        val = cell2mat(val);
    end
end
end
