addpath ../src
addpath ./auxiliary
clear; clc;
%Test the correctness and the history of the error for CRQopt

%%construct matrices
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
if beta == 100
    opts.maxit = 30;
else
    opts.maxit = 130;
end
opts.tol = 1e-15;
opts.method = 2;
opts.checkstep = 1;
opts.returnQ = 1;
%get the solutions of CRQopt by solving QEPmin
[x,infoQEP] = CRQ_Lanczos(A,C,b,opts);

%find the errors and upper bound for the error
[erfQEP,ervQEP] = ErrorHistory(A,infoQEP,vs);

%get the solutions of CRQopt by solving LGopt
opts.method = 1;
[x,infoLGopt] = CRQ_Lanczos(A,C,b,opts);

%find the errors and upper bound for the error
[erfLGopt,ervLGopt] = ErrorHistory(A,infoLGopt,vs);

%plot the convergence history
if beta == 100
   figure(1)
   plot1 = semilogy(ervLGopt,'b-+','linewidth',2); hold on
   plot2 = semilogy(erfLGopt,'k-o','linewidth',2); hold on
   plot3 = semilogy(abs((infoLGopt.mu-lam)/lam),'r-*','linewidth',2); hold on
   
   plot4 = semilogy(ervQEP,'y','linewidth',2); hold on
   plot5 = semilogy(erfQEP,'g','linewidth',2); hold on
   plot6 = semilogy(abs((infoQEP.mu-lam)/lam),'c','linewidth',2); hold on
   
   xlabel('k'); 
   text(4,1e-10,strcat('\beta = ',num2str(beta)),'FontSize',20);
   h = legend([plot2 plot1 plot3 plot5 plot4 plot6],'err_1 by LGopt','err_2 by LGopt','err_3 by LGopt','err_1 by QEPmin','err_2 by QEPmin','err_3 by QEPmin');
   set(h,'FontSize',18);
   set(gca,'FontSize',18)
   axis([0,25,3e-16,5])
   hold off
   print(strcat('correct',num2str(beta),'.eps'),'-depsc')
   
elseif beta == 1000
   figure(1)
   plot1 = semilogy(1:5:size(ervLGopt,2),ervLGopt(1:5:end),'b-+','linewidth',2); hold on
   plot2 = semilogy(1:5:size(erfLGopt,2),erfLGopt(1:5:end),'k-o','linewidth',2); hold on
   plot3 = semilogy(1:5:size(infoLGopt.mu,2),abs((infoLGopt.mu(1:5:end)-lam)/lam),'r-*','linewidth',2); hold on
   
   plot4 = semilogy(1:5:size(ervQEP,2),ervQEP(1:5:end),'y','linewidth',2); hold on
   plot5 = semilogy(1:5:size(erfQEP,2),erfQEP(1:5:end),'g','linewidth',2); hold on
   plot6 = semilogy(1:5:size(infoQEP.mu,2),abs((infoQEP.mu(1:5:end)-lam)/lam),'c','linewidth',2); hold on
   
   xlabel('k'); 
   text(4,1e-10,strcat('\beta = ',num2str(beta)),'FontSize',20);
   h = legend([plot2 plot1 plot3 plot5 plot4 plot6],'err_1 by LGopt','err_2 by LGopt','err_3 by LGopt','err_1 by QEPmin','err_2 by QEPmin','err_3 by QEPmin');
   set(h,'FontSize',18);
   set(gca,'FontSize',18)
   axis([-1,109,3e-16,100])
   hold off
   print(strcat('correct',num2str(beta),'.eps'),'-depsc')
   
else
   figure(1)
   plot1 = semilogy(ervLGopt,'b-+','linewidth',2); hold on
   plot2 = semilogy(erfLGopt,'k-o','linewidth',2); hold on
   plot3 = semilogy(abs((infoLGopt.mu-lam)/lam),'r-*','linewidth',2); hold on
   
   plot4 = semilogy(ervQEP,'y','linewidth',2); hold on
   plot5 = semilogy(erfQEP,'g','linewidth',2); hold on
   plot6 = semilogy(abs((infoQEP.mu-lam)/lam),'c','linewidth',2); hold on
   
   xlabel('k'); 
   text(4,1e-10,strcat('\beta = ',num2str(beta)),'FontSize',20);
   h = legend([plot2 plot1 plot3 plot5 plot4 plot6],'err_1 by LGopt','err_2 by LGopt','err_3 by LGopt','err_1 by QEPmin','err_2 by QEPmin','err_3 by QEPmin');
   set(h,'FontSize',18);
   set(gca,'FontSize',18)
   hold off
end

