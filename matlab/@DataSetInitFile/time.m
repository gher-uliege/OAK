% all time instances where we have observations

function t = time(self)

if ~isempty(self.n)
  t = zeros(1,length(self.n));

  for i=1:length(self.n)
    str = getcopy(self.init,sprintf('Obs%03g.time',self.n(i)));
    t(i) = parseTime(self,str);
  end
else
  i = 1;
  t = [];
  
  while true
    str = get(self.init,sprintf('Obs%03g.time',i),[]);
    
    if isempty(str)
      break;
    else
      t(i) = parseTime(self,str);
      i = i+1;
    end
  end
end
  


function t = parseTime(self,str)
  str = strrep(str,'T',' ');

  %t(i) = mjd(str,31); 
  % yyyy-mm-dd HH:MM:SS  

  [y,m,d,h,mi,s] = datevec(str,31); 
  t = datenum_cal(y,m,d,h,mi,s,self.calendar) - self.torigin;

