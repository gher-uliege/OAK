function [out,coeff] = ndgrid_coeff(X,xi)

tol = 1e-10;

n = length(xi);
M = ones(n+1,n+1);
d = ones(n+1,1);

M(2:n+1,:) = X;
d(2:n+1) = xi';

% find coeff that satisfies
% d = M * coeff

determ = det(M);

if abs(determ) > tol
    % non degenerated case
    coeff = M \ d;
    out = all(0-tol <= coeff & coeff <= 1+tol);
    return;
end

% degenerated case

[U,S,V] = svd(M);

nz = diag(S) > tol;

% remove redundant contraints
U = U(:,nz);
S = S(nz,nz);
V = V(:,nz);

% number of non-redundant constraints
nnz = sum(nz);

% number of coefficients uncontrained
% nuncon + nnz = n+1
nuncon = n+1-nnz;


coeff = zeros(n+1,1);
% is coef initialized with a meanful value?
cinit = false;
best = 0;

z = zeros(nuncon,1);
out = true;

% offset used to convert linear index i to subscribt
% similar to ind2sub in matlab/octave
offset = zeros(nuncon,1);
offset(1) = 1;
for j = 2:nuncon
    offset(j) = offset(j-1) * (n+1);
end

% force one coeff to zero after the other

for i = 1:(n+1)^nuncon
    
    % indices for vertices to test if the corresponding c can be zero
    %[ind(1),ind(2)] = ind2sub((n+1)*ones(nuncon,1),i);
    
    % tmp and ind are here 0-based
    tmp = i - 1;
    
    for j = nuncon:-1:1
        ind(j) = floor(tmp / offset(j));
        tmp = tmp - ind(j) * offset(j);
    end
    % make ind 1-based
    ind = ind+1;
    
    A = zeros(nuncon,n+1);
    
    for j = 1:nuncon
        A(j,ind(j)) = 1;
    end
    
    if any(sum(A,2) == 0)
        % all elements of ind must be different
        % ignore the present case
        continue
    end
    
    %  S * V' * coeff = U' * d
    %  A*c = 0
    
    M2 = [S * V'; A];
    detM2 = det(M2);
    if abs(detM2) < tol
        % still degenerated, look for other possibilities
        continue
    end
    testc = M2 \ [U'* d; z];
    
    %testc
    %all(0 <= testc & testc <= 1)
    
    if all(0-tol <= testc & testc <= 1+tol)
        % yes, this is an interpolation, no extrapolation!
        err = max(abs(M*testc -d));
        
        % but, wait check if it is still a solution to our problem
        % this is necessary if xi is outside a degenerated X, for example
        % xi = [0 1]';
        % X = [0    0.5   1; ...
        %      0      0   0];
        
        
        if err < tol
            % abs(det(M2)) is an indication of the surface of the
            % "triangle", the smaller the better
            if abs(detM2) < best  || ~cinit
                % we have a even better coeff
                coeff = testc;
                cinit = true;
                best = abs(detM2);
                out = false;
            end
        end
    end
end

