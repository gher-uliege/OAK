

tests = {...
   'test_svector',...
   'test_fortranwrite',...
   'test_initfile',...
   'test_model',...
   'test_runoak',...
   'test_runoak2'};

for iindex=1:length(tests);
 
  try
    eval(tests{iindex});
    colordisp('  OK  ','green');
  catch
    colordisp(' FAIL ','red');        
    disp(lasterr)
  end
end


