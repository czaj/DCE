function [LL,g,h] = LL_hmnl_MATlike(Y,Xa,X_str,X_mea,Xmea_exp,err_sliced,EstimOpt,OptimOpt,b0)

LLfun = @(B) LL_hmnl(Y,Xa,X_str,X_mea, Xmea_exp,err_sliced,EstimOpt,B);

if isequal(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        [f,j] = LLfun(b0);
        j(:,EstimOpt.BActive ==0) = 0;
         g = sum(j,1)'; ...
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            h = j'*j;
        end
    else % => EstimOpt.NumGrad == 1 
        f = LLfun(b0);  
        j = numdiff(LLfun,f,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);...
        g = sum(j,1)';   
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            h = j'*j;
        end
    end
else % No gradient
    EstimOpt.NumGrad = 1;
    f = LLfun(b0);   
end

LL = sum(f);