addpath ../src
addpath ./auxiliary
clear;clc;
%Examples where the bound is sharp

%construct matrices
%parameters
n = 1100; m = 100; nm = n-m; zeta = .9;
alpha = 1; beta = 1000; 

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
%true solution
[vs,lam] = CRQ_explicit(A,C,b);

%set the parameters
opts.maxit = 200;
opts.tol = 1e-15;
opts.method = 1;
opts.checkstep = 1;
opts.returnQ = 1;


%get the solutions of CRQopt
[x,info] = CRQ_Lanczos(A,C,b,opts);

%find the errors and upper bound for the error
[erf,erv] = ErrorHistory(A,info,vs);
[boundf,boundv,boundeig] = UpperBound(A,C,lam,info);

%plot the error and the bound
if beta == 100
   figure(1)
   plot1 = semilogy(erv,'b','linewidth',2); hold on
   plot2 = semilogy(boundv,'r-.','linewidth',2);hold on;
   
   xlabel('k'); 
   text(5,1e-10,strcat('\beta = ',num2str(beta)),'FontSize',20);
   h = legend([plot1 plot2],'err_2','Error Bound');
   set(h,'FontSize',16);
   set(gca,'FontSize',18)
   axis([0,31,7e-16,5])
   hold off
   print(strcat('EgArgSharpBeta',num2str(beta),'.eps'),'-depsc')
   
   figure(2)
   plot1 = semilogy(erf,'b','linewidth',2); hold on
   plot2 = semilogy(boundf/(vs'*A*vs),'r-.','linewidth',2);
   xlabel('k'); 
   text(3,1e-10,strcat('\beta = ',num2str(beta)),'FontSize',20);
   h = legend([plot1 plot2],'err_1','Error Bound');
   set(h,'FontSize',16);
   set(gca,'FontSize',18)
   axis([0,16,5e-16,1e3])
   hold off   
   print(strcat('EgValSharpBeta',num2str(beta),'.eps'),'-depsc')
   
   figure(3)
   plot1 = semilogy(abs((info.mu-lam)/lam),'b','linewidth',2); hold on
   plot2 = semilogy(boundeig/abs(lam),'r-.','linewidth',2);
   xlabel('k'); 
   text(3,1e-10,strcat('\beta = ',num2str(beta)),'FontSize',20);
   h = legend([plot1 plot2],'err_3','Error Bound');
   set(h,'FontSize',16);
   set(gca,'FontSize',18)
   axis([0,16,5e-16,1e3])
   hold off   
   print(strcat('EgValSharpBeta',num2str(beta),'.eps'),'-depsc')
   
   
elseif beta == 1000
   figure(1)
   plot1 = semilogy(erv,'b','linewidth',2); hold on
   plot2 = semilogy(boundv,'r-.','linewidth',2);
   xlabel('k'); 
   text(10,1e-10,strcat('\beta = ',num2str(beta)),'FontSize',20);
   h = legend([plot1 plot2],'err_2','Error Bound');
   set(h,'FontSize',16);
   set(gca,'FontSize',18)
   axis([-1,134,7e-16,100])
   hold off
   print(strcat('EgArgSharpBeta',num2str(beta),'.eps'),'-depsc')
   
   figure(2)
   plot1 = semilogy(erf,'b','linewidth',2); hold on
   plot2 = semilogy(boundf/(vs'*A*vs),'r-.','linewidth',2);
   xlabel('k'); 
   text(15,1e-10,strcat('\beta = ',num2str(beta)),'FontSize',20);
   h = legend([plot1 plot2],'err_1','Error Bound');
   set(h,'FontSize',16);
   set(gca,'FontSize',18)
   axis([-1,70,2e-16,1e4])
   hold off   
   print(strcat('EgValSharpBeta',num2str(beta),'.eps'),'-depsc')
   
   figure(3)
   plot1 = semilogy(abs((info.mu-lam)/lam),'b','linewidth',2); hold on
   plot2 = semilogy(boundeig/abs(lam),'r-.','linewidth',2);
   xlabel('k'); 
   text(3,1e-10,strcat('\beta = ',num2str(beta)),'FontSize',20);
   h = legend([plot1 plot2],'err_3','Error Bound');
   set(h,'FontSize',16);
   set(gca,'FontSize',18)
   hold off   
   print(strcat('EgValSharpBeta',num2str(beta),'.eps'),'-depsc')
   
   
else
   figure(1)
   plot1 = semilogy(erv,'b','linewidth',2); hold on
   plot2 = semilogy(boundv,'r-.','linewidth',2);
   xlabel('k'); 
   h = legend([plot1 plot2],'err_2','Error Bound');
   set(h,'FontSize',16);
   set(gca,'FontSize',18)
   hold off
   
   figure(2)
   plot1 = semilogy(erf,'b','linewidth',2); hold on
   plot2 = semilogy(boundf/(vs'*A*vs),'r-.','linewidth',2);
   xlabel('k'); 
   h = legend([plot1 plot2],'err_1','Error Bound');
   set(h,'FontSize',16);
   set(gca,'FontSize',18)
   hold off 
   
   figure(3)
   plot1 = semilogy(abs((info.mu-lam)/lam),'b','linewidth',2); hold on
   plot2 = semilogy(boundeig/abs(lam),'r-.','linewidth',2);
   xlabel('k'); 
   text(3,1e-10,strcat('\beta = ',num2str(beta)),'FontSize',20);
   h = legend([plot1 plot2],'err_3','Error Bound');
   set(h,'FontSize',16);
   set(gca,'FontSize',18)
   %axis([0,16,5e-16,1e3])
   hold off   
   print(strcat('EgValSharpBeta',num2str(beta),'.eps'),'-depsc')
end
