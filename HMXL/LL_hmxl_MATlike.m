function [LL,g,h] = LL_hmxl_MATlike(Y,Xa,Xm,Xs, X_str,X_mea,Xmea_exp, err_sliced,W,EstimOpt,OptimOpt,b0)

% save tmp1
% return

LLfun = @(B) LL_hmxl(Y,Xa,Xm,Xs, X_str,X_mea,Xmea_exp, err_sliced,EstimOpt,B);

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


% j2 = jacobianest(@(B) LL_hmxl(Y,Xa,Xm,Xs, X_str,X_mea,Xmea_exp, err_sliced,EstimOpt,B),b0);