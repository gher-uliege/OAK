function res = oak_loadlog(filename,pat);

fname = [tmpnam '.txt'];
%pat = '^...\.forecast.rms_yo-Hx';
cmd = sprintf('awk ''/%s/ { print $2 } '' "%s" > %s',pat,filename,fname);

%disp(cmd);
[o,ret] = system(cmd);

res = load(fname);

delete(fname);