function [x,ld]=CRQ_explicit(A,C,b)
%%
%----------------------------------------------------------------------
%   Solve CRQ Problem with explicit method
%   
%   Input:
%       A       n-by-n matrix
%       C       n-by-m matrix
%       b       m-dimensional vector
%   Output:
%       x       solution of CRQopt
%

%%
%%
%%Transform the matrices to dence matrices
a = full(A);
c = full(C);
n = min(size(c));

%%QR decomposition and get the maticx PAP
[p,r] = qr(c);
pap = p'*a*p;% the matrix pap
m = size(a,1);
c1 = pap(n+1:m,n+1:m);% the lower right block of PAP

%Eigen decomposition
[q,d,~] = eig(c1);
[delta,pm] = sort(real(diag(d)));
q = real(q(:,pm));% permutation of the eigenvecor, the left to the right is corresponding to the smallest to the largest eigenvalue
gamma = pap(n+1:m,1:n);
r = r(1:n,:);
y = (r')\b;
b = -gamma*y;% the vector of TRS
d = q'*b;% represent b in the basis spanned be q
s = sqrt(1-y'*y);% yhe constrain of length of the vector

% Initial value for the iteration of root finder
id = find(abs(d)>1e-10,1,'first');% the first eigenvalue where di nonzero
ld = delta(id)-abs(d(id))/s;% initial value for ld
ldp = ld+1;% initial value for the ld in the prevoud iteration

%Root finder
itn = 0;
while ldp > ld+1e-10
    itn = itn+1;% iteration number
    ldp = ld;% update ld in the previous iteraion
    ld = ld-2*((f(ld,delta,d,s)+s^2)/fp(ld,delta,d,s))*...
        (sqrt(f(ld,delta,d,s)+s^2)/s-1);% upadate ld
end

%Find x
c = pap(1+n:m,1+n:m);
z = (c-ld*eye(size(c)))\b;
x = [y;z];
x = p*x;
