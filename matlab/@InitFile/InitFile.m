
function retval = InitFile(filename)
keys=cell(1);
values=cell(1);
self.filename = filename;

%args : initifile, key, default ; returned : val
key='Obs.\w*.time';
fid=fopen(filename,'r');
tline=fgetl(fid);
i=1;
while ischar(tline)  %Read all the file and search for the pattern Obs*anything*.time
    tline=fgetl(fid);
    if ~ischar(tline)
        break;
    end
    %val = [];
    %error opening file
      if fid == -1
      error(['unable to open ' filename]);
      end

      line = strtrim(tline);
      
      if isempty(line)
        continue
      elseif line(1) == '!' || line(1) == '#'
        continue
      elseif isempty(regexp(tline,key,'match'));
          continue
      else
        p = strsplit(line,'=');
        
         if length(p) == 2
          keys{i} = strtrim(p{1});
          values{i} = strtrim(p{2});
          
          

            if values{i}(1) == '[' && values{i}(end) == ']'
              values{i}(1) = '{';
              values{i}(end) = '}';
            end
          
            eval(['values{i} = ' values{i} ';']);
            if iscell(values{i})
              if isnumeric(values{i}{1})
                values{i} = cell2mat(values{i});
              end
            end
            i=i+1;
          end          
      end
       

          
end  
         
fclose(fid);
  

self.keys=keys;
self.values=values;
retval = class(self,'InitFile');

%pb: won't take last value but first