function varargout = oak_unpack(ml,x,fillvalue)

if (nargin ==  2)
  fillvalue = NaN;
end

sv = statevector_init(ml.masks{:});

varargout = cell(nargout,1);
[varargout{:}] = statevector_unpack(sv,x,fillvalue);
