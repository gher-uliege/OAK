function retval = DataSetInitFile(filename,n)

self.filename = filename;
self.init = InitFile(filename);
self.n = n;

retval = class(self,'DataSetInitFile');