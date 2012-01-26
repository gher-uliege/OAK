% extracts a column with the syntax SV(:,m)


function self = subsasgn(self,idx,rhs)

%whos rhs
%keyboard
assert(strcmp(idx.type,'()'))
assert(strcmp(idx.subs{1},':'))

m = idx.subs{2};

if strcmp(m,':')
    m = 1:size(self,2);
end

v = cell(self.sv.nvar,1);

if isa(rhs,'SVector')
    for j=1:length(m)
        for i=1:self.sv.nvar
            if 1
                self.variables{i} = rhs.variables{i};
            else
                %disp(['save ' name(self,i,m(j))]);
                filename = gfilename(name(self,i,m(j)));
                filename_rhs = gfilename(name(rhs,i,j));
                
                if ~strcmp(realpath(filename),realpath(filename_rhs))
                    symlink(realpath(filename_rhs),realpath(filename),'delete');
                end
            end
            %disp(['save--- ' name(self,i,m(j)) ' -- ' filename_rhs]);
        end
    end
else
    for j=1:length(m)
        [v{:}] = statevector_unpack(self.sv,rhs(:,j));
        
        for i=1:self.sv.nvar
            filename = name(self,i,m(j));
            %disp(['save ' filename]);
            
            % make parent dir if necessary
            dn = dirname(filename);
            if exist(dn,'dir') ==  0
                syscmd('mkdir -p "%s"',dn);
            end
            
            gwrite(filename,v{i})
        end
    end
end
