function display(self)

sprintf('variables (in %s): \n',self.path)
for i=1:length(self.variables)
  disp(self.variables{i})
end
disp('Members: ')
self.members