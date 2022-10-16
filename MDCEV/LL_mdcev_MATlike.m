function [f,g,h]= LL_mdcev_MATlike(data,EstimOpt,OptimOpt,b0)
% Function calculates loglikelihood of the MDCEV model.
% Return values:
%   f -- loglikelihood value
%   g -- gradient
%   h -- hessian

% save res_LL_mdcev_MATlike
% return

probsfun = @(B) probs_mdcev(data,EstimOpt,B);

if isequal(OptimOpt.GradObj,'on') 
    logprobs = probsfun(b0);
    if isequal(OptimOpt.Hessian,'user-supplied') == 1
        j = -numdiff(probsfun,logprobs,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
        j = j.*data.W;
        g = sum(j,1)';
        h = j'*j;
    else
        j = -numdiff(probsfun,logprobs,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
        j = j.*data.W;
        g = sum(j,1)';
    end
    logprobs = logprobs.*data.W;
    f = -sum(logprobs);
else
    EstimOpt.NumGrad = 1;
    logprobs = probsfun(b0);
    logprobs = logprobs.*data.W;
    f = -sum(logprobs);   
end
