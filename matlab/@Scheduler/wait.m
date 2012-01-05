% wait unil job returns

function wait(self,job)

while isrunning(self,job)
  system('sleep 5');
end