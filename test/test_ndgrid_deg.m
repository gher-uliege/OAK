tol = 1e-7;

% 2d example
X = [0    0.5   1; ...
     0      0   0];

xi = [0 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == false)
assert(max(abs([1; 0; 0] - coeff)) < tol)

xi = [1 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == false)
assert(max(abs([0; 0; 1] - coeff)) < tol)

xi = [0.5 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == false)
%assert(max(abs([.5; 0; .5] - coeff)) < tol)
assert(max(abs([0; 1; 0] - coeff)) < tol)

xi = [0.49 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == false)
assert(max(abs([.02; .98; 0] - coeff)) < tol)

xi = [0.51 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == false)
assert(max(abs([0; .98; 0.02] - coeff)) < tol)


xi = [.2 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == false)
assert(max(abs([0.6; 0.4; 0] - coeff)) < tol)

xi = [2 0]';
[out,~] = ndgrid_coeff(X,xi);
assert(out == true)

xi = [-2 0]';
[out,~] = ndgrid_coeff(X,xi);
assert(out == true)

xi = [0 1]';
[out,~] = ndgrid_coeff(X,xi);
assert(out == true)


% 3d example

% 1-④                   ③
%   |
%   |
%   |         
%   |
%   |
% 0-①-------------------②                
%   |                    |
%   0                    1
%

X = [0  1  1  0; ...
     0  0  1  1; ...
     0  0  0  0];

xi = [0 0 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == false)
assert(max(abs([1; 0; 0; 0] - coeff)) < tol)

xi = [1 0 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == false)
assert(max(abs([0; 1; 0; 0] - coeff)) < tol)

xi = [.5 0 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == false)
assert(max(abs([0.5; 0.5; 0; 0] - coeff)) < tol)

xi = [1 1 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == false)
assert(max(abs([0; 0; 1; 0] - coeff)) < tol)

xi = [.5 .5 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == false)
% 2 solutions
assert(max(abs([0; .5; 0; .5] - coeff)) < tol ||  ...
       max(abs([.5; 0; .5; 0] - coeff)) < tol)

xi = [2 0 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == true)

xi = [0 0 1]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == true)


% 2-|         ③
%   |
%   |
% 1-|         ④
%   |
%   |
% 0-①-------------------②                
%   |         |          |
%   0         1          2
%

% 3d example
X = [0  2  1  1; ...
     0  0  2  1; ...
     0  0  0  0];

xi = [0 0 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == false)
assert(max(abs([1; 0; 0; 0] - coeff)) < tol)

xi = [1 1.5 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == false)
assert(max(abs([0; 0; 0.5; 0.5] - coeff)) < tol)


% 3d example (degenerated in 2 dimensions)
X = [0  1  2  3; ...
     0  0  0  0; ...
     0  0  0  0];

xi = [0 0 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == false)
assert(max(abs([1; 0; 0; 0] - coeff)) < tol)

xi = [.5 0 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == false)
assert(max(abs([0.5; 0.5; 0; 0] - coeff)) < tol)

xi = [2.5 0 0]';
[out,coeff] = ndgrid_coeff(X,xi);
assert(out == false)
assert(max(abs([0; 0; 0.5; 0.5] - coeff)) < tol)



