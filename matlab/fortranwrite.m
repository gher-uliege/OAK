function str = fortranwrite(fmt,varargin)

args = varargin;
iargs = 1;

if fmt(1) == '(' && fmt(end) == ')'
  % fortran format
  
  %part = strsplit(fmt(2:end-1),',');
  part = {};

  i = 2;  
  while i < length(fmt)-1
    [part{end+1},i] = oak_nexttoken(fmt(1:end-1),i,',');
  end
  
  
  for i=1:length(part)
    
    if (part{i}(1) == '"' &&  part{i}(end) == '"') || ...
      (part{i}(1) == '''' &&  part{i}(end) == '''') 
      p{i} = part{i}(2:end-1);
    elseif part{i}(1) == 'I'
      num = strsplit(part{i}(2:end),'.');
      width = num2str(num{1});
      digit = num2str(num{2});
          
      p{i} = sprintf(['%0' num2str(width)  'g'],args{iargs});
      iargs = iargs+1;
      
    else
      error(['I do not know what to do with ' part{i} '. Check ' mfilename ' to see if you can improve it']);
    end
  end
else
  error(['does not look like a fortran format:' fmt]);
end

%  strjoin from octave/matlab
str = strjoin(p,'');
