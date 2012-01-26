% save(self,filenames)
% makes sure that the state vector is available from the filenames

function save(self,path,v,operation)

if nargin == 3
    operation = 'link';
end

if length(v) ~= length(self.variables)
    error('wrong number of elements in v');
end


for j=1:size(self,2)
    for i=1:length(self.variables)
        target = fullfile(path,fmtprint(v{i},j));
        
        if ischar(self.variables{i})
            filename = gfilename(name(self,i,j));
            filenamet = gfilename(target);
            
            if ~strcmp(realpath(filename),realpath(filenamet))
                
                if strcmp(operation,'link')
                    % make a symbolic link (delete target if it exists)
                    symlink(filename,filenamet,'delete');
                elseif strcmp(operation,'copy')
                    copyfile(filename,filenamet);
                else
                    error(['unknown operation ' operation]);
                end
            end
        else
            gwrite(target,self.variables{i});
        end
    end
end
