function [out,coeff] = ndgrid_coeff(X,xi)

tol = 1e-10;

n = length(xi);
M = ones(n+1,n+1);
M2 = ones(n+1,n+1);
d = ones(n+1,1);

% average all X and make vectors relative to this location
% for better numerical precision

xc = mean(X,2);
M(2:n+1,:) = X - repmat(xc,[1 n+1]);
d(2:n+1) = xi' - xc';

% find coeff that satisfies
% d = M * coeff

determ = det(M);

if abs(determ) > tol
    % non degenerated case
    coeff = M \ d;
    out = ~all(0-tol <= coeff & coeff <= 1+tol);
    return;
end

% degenerated case

[U,S,V] = svd(M);

S = diag(S);
nz = find(S > tol);

% number of non-redundant constraints
nnz = length(nz);

% number of coefficients uncontrained
% nuncon + nnz = n+1
nuncon = n+1-nnz;


coeff = zeros(n+1,1);
% is coef initialized with a meanful value?
cinit = false;
best = 0;

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
    
    M2 = zeros(n+1,n+1);
    
    for j = 1:nuncon
        M2(j,ind(j)) = 1;
    end
    
    if any(sum(M2(1:nuncon,:),2) == 0)
        % all elements of ind must be different
        % ignore the present case
        continue
    end
    
    % with S,V
    % % remove redundant contraints
    %  S * V' * coeff = U' * d
    %  A*c = 0
    
    %M2(1:nuncon,:) = A;
    M2(nuncon+1:end,:) = diag(S(nz)) * V(:,nz)';
    
    detM2 = det(M2);
    
    if abs(detM2) < tol
        % still degenerated, look for other possibilities
        continue
    end
    
    d2 = zeros(n+1,1);
    d2(nuncon+1:end) = U(:,nz)'* d;
    
    testc = M2 \ d2;
    
    if ~all(0-tol <= testc & testc <= 1+tol)
        % no, this is a extrapolation!
        % go to next iteration
        continue
    end
        
    err = max(abs(M*testc - d));
    
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

