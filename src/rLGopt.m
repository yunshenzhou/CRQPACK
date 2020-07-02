function [mu,x] = rLGopt(T,w,gamma2,tol)
%
%----------------------------------------------------------------------
%   Purpose
%   Solve the optimization problem 
%   
%   min lam
%    st (T-mu*I)x = -w
%       \|x\| = gamma
%
%----------------------------------------------------------------------
%   
%   Input:
%       T       matrix T
%       w       vector w
%       gamma2  the square of parameter gamma
%       tol     tolerance for iteration of the root finder
%   Output:
%       mu      computed Lagrange multiplier mu
%       x       computed vector x
%
%----------------------------------------------------------------------

k = size(T,1);

%Eigen decomposition
[q,d,~] = eig(T);
[theta,pm] = sort(real(diag(d)));
q = real(q(:,pm));
xi = q'*w;
gamma = sqrt(gamma2);
id = find(abs(d)>1e-10,1,'first');% the first index where xi nonzero

% Initial value for the iteration of root finder
delta0 = sqrt(sum(xi.^2))/gamma;
alpha = theta(id)-delta0;
beta = theta(id);
eta = gamma2-sum(xi(id+1:end).^2./((alpha-theta(id+1:end)).^2));
if eta > 0
    mu = theta(id)-abs(theta(id))/sqrt(eta);
else
    mu = theta(id)-delta0/2;
end


%Root finder
itn = 0;%iteration count
muprevious = mu+1;% initial value for the mu in the prevoud iteration
while abs(muprevious-mu) > tol
    itn = itn+1;% iteration number
    muprevious = mu;% update mu in the previous iteraion
    
    %update the interval (alpha,beta)
    f = sum(xi.^2./((mu-theta).^2))-gamma2;
    if f > 0
        beta = mu;
    else
        alpha = mu;
    end  
    
    %compute values a and b
    a = (mu-theta(id))^3*sum(xi.^2./((mu-theta).^3));
    b = (mu-theta(id))*sum(xi.^2./((mu-theta).^3))-f;
    
    %update the root mu
    if b > 0 
        mu1 = theta(id)-sqrt(a/b);
        if alpha < mu1 && mu1 < beta
            mu = mu1;
        else
            mu = (alpha+beta)/2;
        end
    else
        mu = (alpha+beta)/2;
    end
end

%Find x
x = -(T-mu*eye(k))\w;