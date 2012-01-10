% formatted printing supporting C and Fortran formats

function str = fmtprint(fmt,varargin)

if fmt(1) == '('
  str = fortranwrite(fmt,varargin{:});
else
  str = sprintf(fmt,varargin{:});
end
