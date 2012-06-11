% all time instances where we have observations

function t = time(self)


t = zeros(1,length(self.n));

for i=1:length(self.n)
  str = getcopy(self.init,sprintf('Obs%03g.time',self.n(i)));
  str = strrep(str,'T',' ');

  %t(i) = mjd(str,31); 
  % yyyy-mm-dd HH:MM:SS  
     
  [y,m,d,h,mi,s] = datevec(str,31); 
  t(i) = datenum_cal(y,m,d,h,mi,s,self.calendar) - self.t0;
end
  
