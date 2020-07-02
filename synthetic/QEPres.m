addpath ../src
addpath ./auxiliary
clear; clc;
%Test the residual of the QEP

%construct matrices
%parameters
n = 1100; m = 100; nm = n-m; zeta = .9;
alpha = 1; beta = 100; 

%matrix H and g0
omega = (beta-alpha)/2; xi = -(alpha+beta)/(beta-alpha);
nodes0 = cos((0:nm-1)'*(pi/(nm-1))); nodes = omega*(nodes0-xi);
H = diag(nodes);
theta_min = min(nodes); theta_max = max(nodes);
g0 = ones(nm,1);

%construct A, C and b from H and g0
rng(11);
C = randn(n,m);
[S,R] = qr(C); R = R(1:m,:); 
a = randn(m,1); a = a/(zeta*norm(a));
tb = zeta^2*a; b = R'*tb; A12 = g0*a';
eta = g0'*(g0./nodes)/zeta^2; 
S = [S(:,m+1:n) S(:,1:m)];
A = S*[H,  A12; A12', eta*eye(m)]*S';  
A = 0.5*(A+A');  % adjust to make sure its symmetry

%set the parameters
if beta == 100
    opts.maxit = 30;
else
    opts.maxit = 130;
end
opts.tol = 1e-15;
opts.method = 2;
opts.checkstep = 1;
opts.resopt = 1;
%get the true residuals in each iteraion
[x info] = CRQ_Lanczos(A,C,b,opts);
res = info.res;

%get the residual bound in each iteraion
opts.resopt = 0;
[x info] = CRQ_Lanczos(A,C,b,opts);
resbound = info.res;

%true eigenvalue
[vs lam] = CRQ_explicit(A,C,b);

%plot the residual history
if beta == 100
   figure(1)
   plot1 = semilogy(res,'b','linewidth',2); hold on
   plot2 = semilogy(resbound,'r-.','linewidth',2);
   xlabel('k'); 
   text(5,1e-10,strcat('\beta = ',num2str(beta)),'FontSize',20);
   h = legend([plot1 plot2],'Residual','Residual Bound');
   set(h,'FontSize',18);
   set(gca,'FontSize',18)
   axis([0,31,1e-17,5])
   hold off
   print(strcat('resbound',num2str(beta),'.eps'),'-depsc')
   
   
   
elseif beta == 1000
   figure(1)
   plot1 = semilogy(res,'b','linewidth',2); hold on
   plot2 = semilogy(resbound,'r-.','linewidth',2);
   xlabel('k'); 
   text(10,1e-10,strcat('\beta = ',num2str(beta)),'FontSize',20);
   h = legend([plot1 plot2],'Residual','Residual Bound');
   set(h,'FontSize',18);
   set(gca,'FontSize',18)
   axis([-1,134,1e-17,100])
   hold off
   print(strcat('resbound',num2str(beta),'.eps'),'-depsc')
   
   
else
   figure(1)
   plot1 = semilogy(res,'b'); hold on
   plot2 = semilogy(resbound,'r-.');
   xlabel('k'); 
   h = legend([plot1 plot2],'Residual','Residual Bound');
   set(gca,'FontSize',18)
   print(strcat('resbound',num2str(beta),'.eps'),'-depsc')
   hold off
   
   
end
