function val = gettime(self,key,default)

if nargin == 2
  t = gettime(self.init,key,self.calendar);
else
  t = gettime(self.init,key,self.calendar,default);
end

val = t - self.torigin;
