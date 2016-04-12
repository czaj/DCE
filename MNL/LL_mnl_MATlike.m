function [f,g,h]= LL_mnl_MATlike(Y,Xa,Xm,Xs,W,EstimOpt,OptimOpt,b0)

% save res_LL_mnl_MATlike
% return

LLfun = @(B) LL_mnl(Y,Xa,Xm,Xs,EstimOpt,B);
% [f,j] = LLfun(b0);
%         j(:,EstimOpt.BActive ==0) = 0;
%         g = sum(j,1)'; ...
% EstimOpt.NumGrad = 1;
% LLfun2 = @(B) LL_mnl(Y,Xa,Xm,Xs,EstimOpt,B);
% f2 = LLfun2(b0);  
%         j2 = numdiff(LLfun2,f2,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);...
%         g2 = sum(j2,1)';  
% [g,g2, abs(g-g2)]
% EstimOpt.NumGrad = 0;

%pause;

if isequal(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            [LL, j] = LLfun(b0);
            LL = LL.*W;
            j = j.*W(:, ones(1, size(j,2)));
            f = -sum(LL);
            j(:,EstimOpt.BActive ==0) = 0;
            g = -sum(j); ...
            h = j'*j;
        else
            [LL, j] = LLfun(b0);
            LL = LL.*W;
            j = j.*W(:, ones(1, size(j,2)));
            f = -sum(LL);
            j(:,EstimOpt.BActive ==0) = 0;
            g = -sum(j); ...
        end
    else         
        LL = LLfun(b0);
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            j = -numdiff(LLfun,LL,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);...
            j = j.*W(:, ones(1, size(j,2)));
            g = sum(j,1)'; ...
            h = j'*j;
        else
            j = -numdiff(LLfun,LL,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);...
            j = j.*W(:, ones(1, size(j,2)));
            g = sum(j,1)';
        end
        LL = LL.*W;
        f = -sum(LL);
    end
else
    EstimOpt.NumGrad = 1;
    LL = LLfun(b0);
    LL = LL.*W;
    f = -sum(LL);
    
end