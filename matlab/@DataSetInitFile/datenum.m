function dn = datenum(self,varargin)

dn = datenum_cal(varargin{:},self.calendar) - self.torigin;

