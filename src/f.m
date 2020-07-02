function f = f(ld,delta,d,s)
%
%----------------------------------------------------------------------
%   The function value of the explicit secular equation
%   
%   Input:
%       ld: the independent variable for f
%       delta: a vector of all the eigenvalues of the matrix
%       d: a vector of the numeritors
%       s: the norm of the vector
%   Output:
%       x: the function value f
%

a = d./(delta-ld);%each term
k = find(d==0);% if the numerator of a term is zero, set it to 0
a(k) = 0;
f = sum(a.^2)-s^2;%sum each term and minus s^2
return