function oak_plot_log(assimlog,freelog,var)


rmsfree = oak_loadlog(freelog,['^...\.' var '\.free.rms_yo-Hx']);

rmsf = oak_loadlog(assimlog,['^...\.' var '\.forecast.rms_yo-Hx']);
rmsa = oak_loadlog(assimlog,['^...\.' var '\.analysis.rms_yo-Hx']);
t = oak_loadlog(assimlog,'^...\.temp\.analysis.mjd');
t = t+678942;

rmsfa = [rmsf'; rmsa'];
tfa = [t'; t'];


plot(t,rmsfree,'k-',...
     tfa(1,:),rmsfa(1,:),...
     'ro',tfa(2,:),rmsfa(2,:),'gx',...
     tfa(:),rmsfa(:),'-');

datetick('x');

legend('Free','Forecast','Analysis');