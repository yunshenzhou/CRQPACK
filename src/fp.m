function fp = fp(ld,delta,d,s)
%%
%----------------------------------------------------------------------
%   The derivative of the explicit secular equation
%   
%   Input:
%       ld: the independent variable for f
%       delta: a vector of all the eigenvalues of the matrix
%       d: a vector of the numeritors
%       s: the norm of the vector
%   Output:
%       x: tthe derivative fp
%

a = (d.^2)./((delta-ld).^3);% each term
k = find(d==0);% if the numerator of a term is zero, set it to 0
a(k) = 0;
fp = 2*sum(a);% sum each term to get the derivative
return