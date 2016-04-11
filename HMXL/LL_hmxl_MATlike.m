function [LL,g,h] = LL_hmxl_MATlike(Y,Xa,Xm, X_str,X_mea,Xmea_exp, err_sliced,W,EstimOpt,OptimOpt,b0)

LLfun = @(B) LL_hmxl(Y,Xa,Xm, X_str,X_mea,Xmea_exp, err_sliced,EstimOpt,B);
% tic
% [f,j] = LLfun(b0);
% toc
%         j(:,EstimOpt.BActive ==0) = 0;
%         g = sum(j,1)'; ...
% EstimOpt.NumGrad = 1;
% LLfun2 = @(B) LL_hmxl(Y,Xa,Xm, X_str,X_mea,Xmea_exp, err_sliced,EstimOpt,B);
% tic
% f2 = LLfun2(b0);  
%         j2 = numdiff(LLfun2,f2,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);...
% toc
%         g2 = sum(j2,1)';  
% [g,g2, abs(g-g2)]
% EstimOpt.NumGrad = 0;
%pause;
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