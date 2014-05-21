% accum: 'monthly' for monthly accumulation or function hander
% return index based on time 

function oak_plot_log(assimlog,freelog,var,accum)

if nargin == 3
    accum = 'none';
end

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

% agregation
if ~strcmp(accum,'none')
    
    if strcmp(accum,'monthly')
        dv = datevec(t);
        index = dv(:,2);
    else
        index = accum(t);
    end
    
    rmsfree = accumarray_nanmean(index,rmsfree)';
    rmsf = accumarray_nanmean(index,rmsf)';
    rmsa = accumarray_nanmean(index,rmsa)';
    t = 1:length(rmsf);
end

rmsfa = [rmsf; rmsa];
tfa = [t; t];




%whos rmsf rmsfa
%rg(t)

plot(t,rmsfree,'k-',...
     tfa(1,:),rmsfa(1,:),...
     'ro',tfa(2,:),rmsfa(2,:),'gx',...
     tfa(:),rmsfa(:),'-');

if strcmp(accum,'none')
  datetick('x');
end

legend('Free','Forecast','Analysis');


mfree = sqrt(nanmean(rmsfree.^2));
mf = sqrt(nanmean(rmsf.^2));
ma = sqrt(nanmean(rmsa.^2));

fprintf('%s: mean RMS free %10g forecast %10g analysis %10g \n',var,mfree,mf,ma)


