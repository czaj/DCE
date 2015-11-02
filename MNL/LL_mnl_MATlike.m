function [f,g,h]= LL_mnl_MATlike(Y,Xa,Xs,EstimOpt,OptimOpt,b0)

% save res_LL_mnl_MATlike
% return

LLfun = @(B) LL_mnl(Y,Xa,Xs,EstimOpt,B);

if isequal(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            [LL, j] = LLfun(b0);
            f = -sum(LL);
            j(:,EstimOpt.BActive ==0) = 0;
            g = -sum(j); ...
            h = j'*j;
        else
            [LL, j] = LLfun(b0);
            f = -sum(LL);
            j(:,EstimOpt.BActive ==0) = 0;
            g = -sum(j); ...
        end
    else         
        LL = LLfun(b0);
        f = -sum(LL);
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            j = -numdiff(LLfun,LL,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);...
            g = sum(j,1)'; ...
            h = j'*j;
        else
            j = -numdiff(LLfun,LL,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);...
            g = sum(j,1)';
        end
    end
else
    EstimOpt.NumGrad = 1;
    LL = LLfun(b0);
    f = -sum(LL);
    
end