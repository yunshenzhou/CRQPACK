function [v,info] = CRQ_Lanczos(A,C,b,opts)
%
%----------------------------------------------------------------------
%   Solve CRQ Problem
%
%       min     v'Av
%       st      v'v=1
%               C'v=b
%
%   with projected QEP method.
%   
%   Input:
%       A       n-by-n matrix
%       C       n-by-m matrix
%       b       m-dimensional vector
%       opts    parameters for running
%           opts.maxit: maximum number of Lanczos iteraions 
%               default: min(400,size(A,1))
%           opts.minit: minimum number of Lanczos iteraions 
%               default: 1
%           opts.tol: tolerance of relative residual
%               default: 1e-8
%           opts.checkstep: the number of Lanczos steps between solving 
%               rLGopt or rQEPmin and checking the residual
%               default: 1 
%           opts.return Q: 0 (default): don't return Q_k in structure info 
%                          1: return Q_k in structure info 
%           opts.method: method to solve the problem
%               1 (default): Call Algorithm 2, which uses reduction of 
%                   LGopt and solve rLGopt
%               2: Call Algorithm 3, which uses reduction of QEPmin
%                   and solve rQEPmin    
%           opts.resopt: option for the residual (only valid when method=2)
%               0 (default): using bound to estimate residual
%               1: using exact residual    
%
%   Output:
%       v       solution of CRQ
%       info    information for the computing process
%            info.n0 vector n0
%            info.b0 vector b0
%            info.gamma2 the square of the norm of u gamma^2
%            info.k the number of Lanczos iteration
%            info.T tridiagonal matrix            
%            info.mu Lagrange multipliers (method=1) 
%               or eigenvalues (method=2) in each iteraion
%            info.res upper bound relative resiual of Lagrange equations
%               (method=1) or QEP (method=2) in each iteraion
%            info.Q matrix Q (only avaible when opts.returnQ=1
%
%           if you choose method=2, then additional infomation is returned:
%            info.s the elements of the cell represent eigenvector 
%                   corresponding to the desired eigenvalue of matrix 
%                   [T,-I;-e1*e1'/gamma2,T] in each iteration
%           if you choose method=1, then additional infomation is returned:
%            info.x the elements of the cell represent the solution of 
%                   rLGopt in each iteraions
%
%----------------------------------------------------------------------

%not enough inputs
if nargin < 3
    error(message('MATLAB:min_QEPPL:NotEnoughInputs'));
end

%check the variable opts
if nargin < 4
    maxit = min(400,size(A,1));
    tol = 1e-8;
    minit = 1;
    checkstep = 1;
    method = 1;
else
    if isfield(opts, 'maxit')
        maxit = opts.maxit;
    else
        maxit = min(400,size(A,1));
    end

    if isfield(opts, 'tol')
        tol = opts.tol;
    else
        tol = 1e-8;
    end

    if isfield(opts, 'minit')
        minit = opts.minit;
    else
        minit = 1;
    end

    if isfield(opts, 'checkstep')
        checkstep = opts.checkstep;
    else
        checkstep = 1;
    end
    if isfield(opts, 'method')
        method = opts.method;
    else
        method = 1;
    end
    if isfield(opts, 'resopt')
        resopt = opts.resopt;
    else
        resopt = 0;
    end
    if isfield(opts, 'returnQ')
        returnQ = opts.returnQ;
    else
        returnQ = 0;
    end
end

%find vector n0
n0 = C*((C'*C)\b);
gamma2 = 1-norm(n0)^2;

%check the existence of solution
zero = 1e-10;
if gamma2 < -zero %no solution
    error('no solution');
elseif gamma2 < zero %there exists a unique solution for constraints
    v = n0;
else %multiple solutions for constraints
    %compute b0
    x = A*n0;
    b0 = x-C*((C'*C)\(C'*x));
    
    %%find solution of QEP
    %construct the function of matrix-vector multiplication PAPx
    Ct = C';
    [L, U, P]  =  lu(full(Ct*C));
    afun = @(x) A*x;
    acfun = @(y) y-C*(U\(L\(P*(Ct*y))));
    fun = @(x) acfun(afun(x));
    
    if method ==  2 %solve QEPmin
        [~,Q,info] = QEPmin(fun,b0,gamma2,norm(A,1),maxit,tol,minit,checkstep,resopt);
        %compute the solution of CRQopt
        s = info.s{info.k};
        y = s(info.k+1:2*info.k);
        w = s(1:info.k);
        e1 = zeros(info.k,1);
        e1(1) = norm(b0);
        v = -Q*(gamma2/(norm(b0)*w(1)))*y+n0;
    else %solve LGopt
        [x,Q,info] = LGopt(fun,b0,gamma2,norm(A,1),maxit,tol,minit,checkstep);        
        %compute the solution of CRQopt
        v = Q*x+n0;
    end
end

%Return the structure info
info.b0 = b0;
info.n0 = n0;
info.gamma2 = gamma2;
if returnQ
    info.Q = Q;
end
