% full(self)
% load SVector in memory

function x = full(self)

v = {};

for i=1:length(self.variables)  
  if ischar(self.variables{i})
    
    if isempty(self.members) 
      % load single state
      v{i} = gread(fullfile(self.path,self.variables{i}));
    else
      % load ensemble of states
      for j=1:size(self,2)
        %disp(['load ' name(self,i,j)]);
        tmp{j} = gread(name(self,i,j));
        %size(tmp{i})
        
        %if j==self.members(1)
          % make allocation
        %  v{i} = zeros(numel(tmp),length(self.members));
        %end
        %whos tmp
        %size(v{i})
        
        %v{i}(:,j) = tmp(:);
      end
      
      v{i} = cat(my_ndims(tmp{j})+1,tmp{:});
    end    
  else    
    % copy state which is already in memory
    v{i} = self.variables{i};
  end  
end

x = statevector_pack(self.sv,v{:});



function d = my_ndims(v)
  if isvector(v)
    d = 1;
  else
    d = ndims(v);
  end
  
