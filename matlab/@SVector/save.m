% save(self,filenames)
% makes sure that the state vector is available from the filenames

function save(self,path,v,operation)

if nargin == 3
    operation = 'link';
end

if length(v) ~= length(self.variables)
    error('wrong number of elements in v');
end


source = {};
dest = {};

for j=1:size(self,2)
    for i=1:length(self.variables)
        target = fullfile(path,fmtprint(v{i},j));
        
        if ischar(self.variables{i})
            filename = gfilename(name(self,i,j));
            filenamet = gfilename(target);
            
            if ~strcmp(realpath(filename),realpath(filenamet))                
                source{end+1} = realpath(filename);
                dest{end+1} = filenamet;
            end
        else
            gwrite(target,self.variables{i});
        end
    end
end

% avoid to copy the same file multiple times
[dest,ii] = unique(dest);
source = source(ii);

for i=1:length(source)       
    if strcmp(operation,'link')
        % make a symbolic link (delete target if it exists)
        symlink(source{i},dest{i},'delete');
    elseif strcmp(operation,'copy')
        fprintf('copy %s to %s\n',source{i},dest{i});
        copyfile(source{i},dest{i});
    else
        error(['unknown operation ' operation]);
    end    
end
