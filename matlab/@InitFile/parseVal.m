function val = parseVal(self,str)

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
