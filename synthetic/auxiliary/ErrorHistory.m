function [erf,erv] = ErrorHistory(A,info,vs)
%
%----------------------------------------------------------------------
%   Purpose
%   Compute the error of the objective function and solution vector
%   in each iteration
%   
%----------------------------------------------------------------------
%   
%   Input:
%       A       n-by-n matrix
%       info    info structure output from CRQopt function
%       vs      the true solution of the CRQ problem
%
%   Output:
%       erf     error of the objective function in each iteration
%       erv     error of the vector in each iteration
%
%----------------------------------------------------------------------

%the true objective function value
MinValue = vs'*A*vs;
if isfield(info,'s') %solved QEPmin
    for k = 1:info.k %for every iteration
        %construct the solution of CRQopt    
        s = info.s{k};
        if size(s,2)>0 %s is nonempty
            y = s(k+1:2*k);
            w = s(1:k);
            u = -info.Q(:,1:k)*(info.gamma2/(norm(info.b0)*w(1)))*y;

            %compute the objective function value
            v = u+info.n0;
            value = v'*A*v;

            %compute the error of objective function value
            erf(k) = abs(value-MinValue)/abs(MinValue);
            erv(k) = norm(vs-v);
        else
            %no information in that iterati%/abs(MinValue);
            erf(k) = NaN;
            erv(k) = NaN;
        end
    end
else %solved LGopt
    for k = 1:info.k %for every iteration
        %construct the solution of CRQopt    
        x = info.x{k};
        if size(x,2) > 0 %s is nonempty
            v = info.Q(:,1:k)*info.x{k}+info.n0;

            %compute the objective function value
            value = v'*A*v;

            %compute the error of objective function value
            erf(k) = abs(value-MinValue)/abs(MinValue);
            erv(k) = norm(vs-v);
        else
            %no information in that iteration;
            erf(k) = NaN;
            erv(k) = NaN;
        end
    end
end