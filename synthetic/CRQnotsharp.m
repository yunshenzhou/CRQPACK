addpath ../src
addpath ./auxiliary
clear;clc;
%Examples where the bound is not sharp

%construct matrices
%parameters
n = 1100; m = 100; nm = n-m; zeta = .9;
alpha = 2; beta = 1000; 

%matrix H and g0
omega = (beta-alpha)/2; xi = -(alpha+beta)/(beta-alpha);
nodes0 = cos((0:nm-2)'*(pi/(nm-2))); nodes = omega*(nodes0-xi);
nodes = [nodes;1];
H = diag(nodes);
theta_min = min(nodes); theta_max = max(nodes);
g0=exp(-5e-3*(1:nm))';

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


%true solution
[vs,lam] = CRQ_explicit(A,C,b);

%set the parameters
opts.maxit = 200;
opts.tol = 1e-8;
opts.method = 2;
opts.checkstep = 1;
opts.returnQ = 1;
%get the solutions of CRQopt
[x,info] = CRQ_Lanczos(A,C,b,opts);

%find the errors and upper bound for the error
[erf,erv] = ErrorHistory(A,info,vs);
[boundf,boundv,boundeig] = UpperBound(A,C,lam,info);
[boundf2,boundv2,boundeig2] = UpperBound2(A,C,lam,info);

%plot the results
x_ax=[1:5:opts.maxit opts.maxit];
figure(1)
plot1 = semilogy(erv,'b','linewidth',2); hold on
plot2 = semilogy(boundv,'r-.','linewidth',2);
plot3 = semilogy(x_ax,boundv2(x_ax),'k-o','linewidth',2);

xlabel('k'); 
h = legend([plot1 plot2 plot3],'err_2','Error bound by \kappa','Error bound by \kappa_+');
set(h,'FontSize',16,'Location','southwest');
set(gca,'FontSize',18)
axis([0,200,4e-6,1e7])
hold off
print('boundf2.eps','-depsc')

figure(2)
plot1 = semilogy(erf,'b-','linewidth',2); hold on
plot2 = semilogy(boundf/abs(vs'*A*vs),'r-.','linewidth',2);
plot3 = semilogy(x_ax,boundf2(x_ax)/abs(vs'*A*vs),'k-o','linewidth',2);

xlabel('k'); 
h = legend([plot1 plot2 plot3],'err_1','Error bound by \kappa','Error bound by \kappa_+');
set(h,'FontSize',16,'Location','southwest');
set(gca,'FontSize',18)
axis([0,200,1e-8,1e10])
hold off   
print('boundv2.eps','-depsc')
   
figure(3)
plot1 = semilogy(abs((info.mu-lam)/lam),'b','linewidth',2); hold on
plot2 = semilogy(boundeig/abs(lam),'r-.','linewidth',2);
plot3 = semilogy(x_ax,boundeig2(x_ax)/abs(lam),'k-+','linewidth',2);

xlabel('k'); 
h = legend([plot1 plot2 plot3],'err_3','Error bound by \kappa','Error bound by \kappa_+');
set(h,'FontSize',16,'Location','southwest');
set(gca,'FontSize',18)
hold off   
