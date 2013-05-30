% extract next token from string str starting a i
% all tokens are separated by sep
% but takens include always opening and closing quotes and parenthesis.

function [token,i] = oak_nexttoken(str,i,sep);

if nargin == 2
  sep = ',';
end

j = i;

while str(j) ~= sep

  if str(j) == '"'
    j = j+1;   
    while str(j) ~= '"'
      j = j+1;
    end
  end

  if str(j) == '('
    j = j+1;   
    while str(j) ~= ')'
      j = j+1;
    end
  end
  
  if j == length(str)
    break;
  end
  
  j = j+1;
end

if str(j) == sep  
  token = str(i:j-1);
else
  token = str(i:j);
end

i = j+1;
