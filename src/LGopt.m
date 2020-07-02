function [x,Q,info] = LGopt(H,g,gamma2,nmH,maxit,tol,minit,checkstep)
%
%----------------------------------------------------------------------
%   Purpose
%   Solve the Lagrange equations  
%   
%   min lam
%    st (H-lam*I)x=-g
%       \|x\|=gamma
%
%   by Krylov subspace method
%----------------------------------------------------------------------
%   
%   Input:
%       H           matrix H, or a handle for function of multiplication Hx
%       g           vector g
%       gamma2      the square of parameter gamma
%       nmH         norm of matrix H
%       maxit       max Lanczos steps
%       tol         tolerance for Lanczos iteration relative residual<tol
%       minit       min Lanczos steps
%       checkstep   the number of Lanczos steps between solving rLGopt and
%                   check the residual
%   Output:
%       z       computed solution of LGopt
%       info    information for Lanczos iteration
%
%            info.k the number of Lanczos iteration
%            info.T tridiagonal matrix
%            info.Q matrix Q
%            info.mu Lagrange multipliers in each iteraion
%            info.res relative residual in each iteraion
%
%----------------------------------------------------------------------


%Initialize
n = size(g,1);
zero = 1e-6;
Q = zeros(n,maxit+1);
alpha = zeros(maxit,1);
beta = zeros(maxit+1,1);
beta(1) = norm(g);
if beta(1) < zero
    disp('input vector is zero');    
    return
end

%Lanczos iterations
q0 = zeros(n,1);
qk = g/beta(1);
Q(:,1) = qk;
for k = 1:maxit 
    %multiply by H
    if isnumeric(H)
        p = H*qk;
    else
        p = H(qk);
    end

    %orthogonalization
    r = p-beta(k)*q0;
    alpha(k) = qk'*p;
    r = r-alpha(k)*qk;
    beta(k+1) = norm(r);
    
    %check the condition to solve rLGopt
    if (k >= minit && mod(k-minit,checkstep) == 0) || k == maxit
        
        %%solve rLGopt
        %construct the tridiagonal matrix
        T = diag(alpha(1:k))+diag(beta(2:k),1)+diag(beta(2:k),-1);
        e1 = zeros(k,1);
        e1(1) = norm(g);

        %solve rLGopt by solving the secular equations
        [mu(k) x] = rLGopt(T,e1,gamma2,tol);
        xx{k} = x;
        %check the relative residual
        res(k) = abs(beta(k+1)*x(k))/((nmH+abs(mu(k)))*norm(x)+norm(g));
        
        %stopping condition
        if  res(k) < tol || k == maxit || beta(k+1) < zero           
            if k < maxit
                disp(['Projected LGopt Converges']);
            end
            disp(['Projected LGopt runs ' num2str(k) ' iterations.']);     
            %generate the structure info
            Q = Q(:,1:k);
            info.k = k;
            info.T = T;
            info.mu = mu;
            info.gamma2 = gamma2;   
            info.res = res;
            info.x = xx;
            return
        end
    end
    
    %update the basis
    q0 = qk;
    qk = r/beta(k+1);
    Q(:,k+1) = qk;
end
