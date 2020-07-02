function [vbar2,Eigenvectors] = ncut(W,nbEigenValues,B,c, dataNcut)
% function [Eigenvectors,Eigenvalues] = ncut(W,nbEigenValues,dataNcut);
% 
% Input:
%     W= symmetric similarity matrix
%     nbEigenValues=  number of Ncut eigenvectors computed
%     dataNcut= optional parameters
%
%     default parameters for dataNcut:
%     dataNcut.offset = 5e-1; offset in the diagonal of W
%     dataNcut.verbose = 0; 0 for verbose off mode, 1,2,3 for verbose on modes
%     dataNcut.maxiterations = 100; max number of iterations in eigensolver
%     dataNcut.eigsErrorTolerance = 1e-6; error tolerance in eigensolver
%     dataNcut.valeurMin=1e-6; % truncates any values in W less than valeurMin
% 
% Output: 
%    Eigenvectors= continuouse Ncut eigenvectos, size = length(W) x nbEigenValues
%    Eigenvalues= Ncut eigenvalues, size = 1x nbEigenValues
%
% Timothee Cour, Stella Yu, Jianbo Shi, 2004.
if nargin < 2
    nbEigenValues = 8;
end
if nargin < 5
    dataNcut.offset = 5e-1;
    dataNcut.verbose = 0;
    dataNcut.maxiterations = 300;
    dataNcut.eigsErrorTolerance = 1e-8;
    dataNcut.valeurMin=1e-6;
end

% make W matrix sparse
W = sparsifyc(W,dataNcut.valeurMin);
% check for matrix symmetry
if max(max(abs(W-W'))) > 1e-10 %voir (-12) 
    error('W not symmetric');
end

n = size(W,1);
nbEigenValues = min(nbEigenValues,n);
offset = dataNcut.offset;


% degrees and regularization
d = sum(abs(W),2);
dr = 0.5 * (d - sum(W,2));
dr = dr + offset;
Dinvsqrt = 1./sqrt(d+eps);
P = spmtimesd(spdiags(d,0,n,n)-W,Dinvsqrt,Dinvsqrt);
clear W;

options.issym = 1;
     
if dataNcut.verbose
    options.disp = 3; 
else
    options.disp = 0; 
end
options.maxit = dataNcut.maxiterations;
options.tol = dataNcut.eigsErrorTolerance;

options.v0 = ones(size(P,1),1);
options.p = max(35,2*nbEigenValues); %voir
options.p = min(options.p,n);



%linear constraints d^(1/2)v'1 = 0
tp = size(B,1)+1;
B(tp,:) = sqrt(d(:));
c(tp) = 0;
%linear constraints c+ = sqrt(1/b) c- = -sqrt(b)
idx = find(c>0);
idx2 = find(c<0);
B1 = sum(B(idx,:),1);
id1 = find(B1>0);
B1 = sum(B(idx2,:),1);
id2 = find(B1>0);
b = sum(d(id1))/sum(d(id2));
c(idx) = sqrt(1/b)*c(idx);
c(idx2) = sqrt(b)*c(idx2);
%linear constraints d^(1/2)c
for j = 1:tp-1
    id = find(B(j,:)>0);
    c(j) = c(j)*sqrt(d(id));
end

%normalize c
c = c/sqrt(sum(d));

%%% Lanczos method for QEP
%set parameters
opts.maxit = 300;
opts.minit = 120;
opts.tol = 8e-5;
opts.method = 2;
opts.checkstep = 5;
%solve the optimzation problem
[vbar2,info] = CRQ_Lanczos(P,B',c,opts);


 %Discrete Eigenvectors
Eigenvectors = zeros(n,2);
k1 = find(vbar2>0);
k2 = find(vbar2<= 0);
Eigenvectors(k1,1) = 1;
Eigenvectors(k2,2) = 1;

