function [boundf,boundv,boundeig] = UpperBound(A,C,lam,info)
%
%----------------------------------------------------------------------
%   Purpose
%   Compute the error bound of the objective function and solution vector \
%   in each iteration
%   
%----------------------------------------------------------------------
%   
%   Input:
%       A       n-by-n matrix
%       C       n-by-m matrix
%       lam     optimal value for LGopt
%       info    info structure output from CRQopt function
%
%   Output:
%       boundf  error bound of the objective function in each iteration (by
%       kappa)
%       boundv  error bound of the vector in each iteration (by kappa)
%       boundf_kappap error bound of the objective function in each iteration 
%       (by kappa+)
%       boundv_kappap error bound of the vector in each iteration (by kappa+)
%
%----------------------------------------------------------------------


%projection matrix
P = eye(max(size(A)))-C*inv(C'*C)*C';

%eigenvalues of PAP
PAP = P*A*P;
eigPAP = sort(real(eig(PAP)));

%remove some zeros to form eigenvalues of H
[~,pm] = sort(abs(eigPAP),'ascend');
eigPAP(pm(1:size(C,2))) = [];

%select the smallest, second smallest and the largest eigenvalue
tn = eigPAP(end);
t1 = eigPAP(1);
%%compute the upper bounds
hopt = tn-lam;
%error bound in by kappa
kappa = (tn-lam)/(t1-lam);
gammat = (sqrt(kappa)+1)/(sqrt(kappa)-1);
%function value of Chebyshev polynomials
chev = 2./(gammat.^(1:info.k)+gammat.^(-(1:info.k)));

%get the bounds
boundf = 4*hopt*info.gamma2*chev.^2;
boundv = 2*sqrt(kappa)*sqrt(info.gamma2)*chev;
boundeig = (boundf+norm(info.b0)*boundv)/info.gamma2;