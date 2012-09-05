function oak_plot_log(assimlog,freelog,var)


rmsfree = oak_loadlog(freelog,var,'free','rms_yo-Hx');
rmsf = oak_loadlog(assimlog,var,'forecast','rms_yo-Hx');
rmsa = oak_loadlog(assimlog,var,'analysis','rms_yo-Hx');
t = oak_loadlog(assimlog,var,'analysis','mjd');

t = t+678942;

if length(rmsf) ~= length(rmsfree)
  warning('timming time series');
  
  i = min([length(rmsf) length(rmsfree)]);
  rmsf = rmsf(1:i);
  rmsa = rmsa(1:i);
  rmsfree = rmsfree(1:i);
  t = t(1:i);
end

rmsfree = rmsfree(~isnan(t));
rmsf = rmsf(~isnan(t));
rmsa = rmsa(~isnan(t));
t = t(~isnan(t));


rmsfa = [rmsf; rmsa];
tfa = [t; t];

whos rmsf rmsfa

rg(t)
plot(t,rmsfree,'k-',...
     tfa(1,:),rmsfa(1,:),...
     'ro',tfa(2,:),rmsfa(2,:),'gx',...
     tfa(:),rmsfa(:),'-');

datetick('x');

legend('Free','Forecast','Analysis');