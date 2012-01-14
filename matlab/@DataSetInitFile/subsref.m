function val = subsref(self,idx)

assert(strcmp(idx.type,'.'));

if strcmp(idx.subs,'filename')
  val = self.filename;
end
