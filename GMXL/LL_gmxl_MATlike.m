function [LL,g,h]= LL_gmxl_MATlike(YY,XXa,XXm,XXs,XXt,err_sliced,W,EstimOpt,OptimOpt,b0)

% save res_LL_mxl_MATlike;
% return

LLfun = @(B) LL_gmxl(YY,XXa,XXm,XXs,XXt,err_sliced,EstimOpt,B);

if isequal(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        [f,j] = LLfun(b0);
        j(:,EstimOpt.BActive ==0) = 0;
        j = j.*W(:, ones(1,size(j,2)));
        g = sum(j,1)'; ...
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            h = j'*j;
        end
    else % => EstimOpt.NumGrad == 1 
        f = LLfun(b0);  
        j = numdiff(LLfun,f,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);...
        j = j.*W(:, ones(1,size(j,2)));
        g = sum(j,1)';   
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            h = j'*j;
        end
    end
else % No gradient
    EstimOpt.NumGrad = 1;
    f = LLfun(b0);   
end
f = f.*W;
LL = sum(f);
