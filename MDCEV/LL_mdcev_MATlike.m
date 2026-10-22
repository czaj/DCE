function [f,g,h]= LL_mdcev_MATlike(data,EstimOpt,OptimOpt,b0)
% Function calculates loglikelihood of the MDCEV model.
% Return values:
%   f -- loglikelihood value
%   g -- gradient
%   h -- hessian

% save res_LL_mdcev_MATlike
% return

probsfun = @(B) probs_mdcev(data,EstimOpt,B);
W_task = data.W(:);

if isequal(OptimOpt.GradObj,'on') 
    if EstimOpt.NumGrad == 0
        [nll,j] = probsfun(b0);
        j(:,EstimOpt.BActive == 0) = 0;
    else
        nll = probsfun(b0);
        j = numdiff(probsfun,nll,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
    end
    j = j.*W_task; g = sum(j,1)';
    if isequal(OptimOpt.Hessian,'user-supplied'), h = j'*j; end
    f = sum(nll.*W_task);
else
    EstimOpt.NumGrad = 1;
    f = sum(probsfun(b0).*W_task);
end
