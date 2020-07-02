function [z,Q,info] = QEPmin(H,g,gamma2,nmH,maxit,tol,minit,checkstep,resopt)
%
%----------------------------------------------------------------------
%   Purpose
%   Solve the smallest real eigenvalue and corresponding eigenvector of
%   QEP (H-lam I)^2z=gg'/gamma2*z 
%   
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
%       checkstep   the number of Lanczos steps between solving rQEPmin and
%                   checking the residual
%   Output:
%       z       computed eigenvector
%       info    information for Lanczos iteration
%
%            info.k the number of Lanczos iteration
%            info.s the elements of the cell represent eigenvector 
%                   corresponding to the desired eigenvalue of matrix 
%                   [T,-I;-e1*e1'/gamma2,T] in each iteration
%            info.T tridiagonal matrix
%            info.Q matrix Q
%            info.mu eigenvalues in each iteraion
%            info.res relative residual in each iteraion
%
%----------------------------------------------------------------------


%Initialize
n = size(g,1);
zero = 1e-6;
nmg = norm(g);

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
    
    %update the basis
    if beta(k+1) >= zero    
        q0 = qk;
        qk = r/beta(k+1);
        Q(:,k+1) = qk;
    end
    %check the condition to solve rQEPmin
    if k == maxit || (k >= minit && mod(k-minit,checkstep) == 0)
        %%solve rQEPmin
        %construct the linearized matrix for QEP
        T = diag(alpha(1:k))+diag(beta(2:k),1)+diag(beta(2:k),-1);
        e1 = zeros(k,1);
        e1(1) = norm(g);    
        mat = [T,-eye(k);-e1*e1'/gamma2,T];

        %find the minimum eigenvalue
        [v e] = eig(mat);
        e = diag(e);
        [lam(k), id] = min(real(e));

        %compute the eigenvectors
        v = v(:,id);
        s{k} = v;
        %sepearte the eigenvector into two parts
        y = v(k+1:2*k);
        w = v(1:k);    

        %check the relative residual
        if resopt
            %true residual
            if isnumeric(H)
                temp = beta(k+1)*y(k)*qk+beta(k+1)*w(k)*(H*qk-lam(k)*qk);
            else
                temp = beta(k+1)*y(k)*qk+beta(k+1)*w(k)*(H(qk)-lam(k)*qk);
            end
            res(k) = norm(temp)/(((nmH+abs(lam(k)))^2+nmg^2/gamma2)*norm(w));
        else
        %residual bound
        res(k) = (abs(beta(k+1))*(abs(y(k))+(nmH+abs(lam(k)))*abs(w(k))))...
            /(((nmH+abs(lam(k)))^2+nmg^2/gamma2)*norm(w));
        end
        %check the stopping condition
        if  res(k) < tol || k == maxit || beta(k+1) < zero  
            if k < maxit
                disp(['QEPmin Converges']);
            end
            disp(['QEPmin runs ' num2str(k) ' Lanczos iterations.']);
            %generate the computed eigenvector
            z = Q(:,1:k)*w;
            Q = Q(:,1:k);
            %generate the structure info
            info.k = k;
            info.s = s;
            info.T = T;
            info.mu = lam;
            info.gamma2 = gamma2;   
            info.res = res;
            return
        end
    end   
end
