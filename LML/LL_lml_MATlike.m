function [f,g,h]= LL_lml_MATlike(GridProbs,b_mtx,W,EstimOpt,OptimOpt,b0)

% save res_LL_lml_MATlike;
% return

LLfun = @(B) LL_lml(GridProbs,b_mtx,EstimOpt,B);

if isequal(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        [f,j] = LLfun(b0);
        j(:,EstimOpt.BActive == 0) = 0;
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
            j = numdiff(LLfun,f,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
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

end