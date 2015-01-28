% b = dirname(path)
% strip directory from filenames

function d = dirname(path)

if path(end) == filesep
    d = path;
else
    i = find(path == filesep,1,'last');
    
    if isempty(i)
        % no /
        d = '.';
    else
        d = path(1:i-1);
    end
end    
    
