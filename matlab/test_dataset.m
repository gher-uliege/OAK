initfile = 'test_assim.init';
clear time

data = DataSetInitFile(initfile,1:2);
time(data)

assert(strcmp(data.filename,initfile))