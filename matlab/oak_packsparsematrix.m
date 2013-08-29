function [H,valid1,valid2] = oak_packSparseMatrix(Hindex,Hcoeff,ML1,ML2)

valid1 = ones(ML1.effsize,1);
valid2 = ones(ML2.effsize,1);

N = length(Hcoeff);

  % allocation if all points where valid
  Hi = zeros(N,1);
  Hj = zeros(N,1);
  Hs = zeros(N,1);
  val1 = zeros(N,1);
  val2 = zeros(N,1);
   
  for i=1:N
    % space 1: destination
    % transform [Hindex(1,i) Hindex(2,i) Hindex(3,i) Hindex(4,i)] into the
    % linear index linindex1 and trapp error in variable val1

    [Hi(i),val1(i)] = oak_sub2ind(ML1,Hindex(1,i),Hindex(2,i),Hindex(3,i),Hindex(4,i));

    % space 2: origin
    % transform [Hindex(5,i) Hindex(6,i) Hindex(7,i) Hindex(8,i)] into the
    % linear index linindex2 and trapp error in variable val2

    [Hj(i),val2(i)] = oak_sub2ind(ML2,Hindex(5,i),Hindex(6,i),Hindex(7,i),Hindex(8,i));
    
    if val1(i)
      valid1(Hi(i)) = val2(i);      
    end
    
    if val2(i) 
      valid2(Hj(i)) = val1(i);
    end
  end

val = val1 & val2;

rg(Hi(val))
rg(Hj(val))
ML2.effsize,ML1.effsize

H = sparse(Hi(val),Hj(val),Hcoeff(val),ML1.effsize,ML2.effsize);


