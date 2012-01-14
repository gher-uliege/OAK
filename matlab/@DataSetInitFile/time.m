% all time instances where we have observations

function t = time(self)

for i=self.n
  str = get(self.init,sprintf('Obs%03g.time',i));
  str = strrep(str,'T',' ');

  t(i) = mjd(str,31); % yyyy-mm-dd HH:MM:SS  
end
  
