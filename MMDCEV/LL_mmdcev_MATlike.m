function [f,g,h]= LL_mmdcev_MATlike(data,EstimOpt,OptimOpt,b0)
% Function calculates loglikelihood of the MDCEV model.
% Return values:
%   f -- loglikelihood value
%   g -- gradient
%   h -- hessian

% save res_LL_mnl_MATlike
% return

LLfun = @(B) probs_mmdcev(data,EstimOpt,B);

W = reshape(data.W,EstimOpt.NCT,[])';
W = W(:,1);

if isequal(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        [f,j] = LLfun(b0);
        j(:,EstimOpt.BActive == 0) = 0;
%          j2 = numdiff(LLfun,f,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
%          [sum(j)', sum(j2)', abs(sum(j) -sum(j2))']
%         pause
        j = j.*W;
        g = sum(j);
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            h = j'*j;
        end
        
    else % => EstimOpt.NumGrad == 1
        f = LLfun(b0);
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            j = numdiff(LLfun,f,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
            j = j.*W;
            g = sum(j,1)';
            h = j'*j;
        else % OptimOpt.Hessian ~= 'user-supplied'
            j = numdiff(LLfun,f,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);...
            j = j.*W;
            g = sum(j,1)';
        end
    end
else % No gradient
    EstimOpt.NumGrad = 1;
    f = LLfun(b0);
end
f = f.*W;
f = sum(f);