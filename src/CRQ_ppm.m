function [v, erf,erv] = CRQ_ppm(A,C,b,alpha,tolPPM,kPPM,history,vs)
%%
%----------------------------------------------------------------------
%   Purpose
%   Compute the CRQ problem by projected power method
%   
%----------------------------------------------------------------------
%   
%   Input:
%   A: n by n sparse symmetric matrix
%   C: n by m sparse symmetric matrix with rank(C) = m
%   b: n-dimensional vector
%   alpha: make H = alpha*speye(n)-A pos. semi-def. 
%   tolPPM: PPM tolerance (not used, reserved for test)
%   kPPM: PPM max iterations
%   history: 1/0 indicator of whether to store the \|v_k-v\| and
%   |v_k'Av_k-v'Av| for every iteraion or not
%   vs: if you assign history=1, then you need to provide the true solution
%   of CRQ problem v
%
%   Output:
%   v: solution of CRQ
%   erf: \|v_k-v\| for each k if you choose history=1
%   erv: |v_k'Av_k-v'Av| for each k if you choose history=1
%
%----------------------------------------------------------------------


n0 = C*((C'*C)\b);

gamma = sqrt(1-norm(n0)^2);

v = n0;
Ct = C';
InvCtC = Ct*C;

for k = 1:kPPM

    v1 = v; 
    %multiply by A
    w = alpha*v-A*v;
    %multiply by P
    y = InvCtC\(Ct*w);
    x = w-C*y;
    %normalize
    v = gamma*x/norm(x) + n0;
    
    %record the history of error
    if history == 1
        erf(k) = norm(v-vs);
        erv(k) = abs(v'*A*v-vs'*A*vs);
    end
        
    %stopping condition
    if norm(v - v1) < tolPPM
       disp('PPM converges to tolPPM');
       break;
    end

end
disp(['PPM runs ' num2str(k) ' iterations.']);
if history == 0
    erf = [];
    erv = [];
end

end
