function [f,g,h]= LL_gwmnl_MATlike(Y,Xa,Weights,EstimOpt,OptimOpt,b0)

LLfun = @(B) LL_gwmnl(Y,Xa,EstimOpt,B);
% [f,j] = LLfun(b0);
%         g = sum(j,1)'; ...
% EstimOpt.NumGrad = 1;
% LLfun2 = @(B) LL_gwmnl(Y,Xa,EstimOpt,B);
% f2 = LLfun2(b0);  
%         j2 = numdiff(LLfun,f2,b0,isequal(OptimOpt.FinDiffType,'central'),ones(EstimOpt.NVarA,1));...
%         g2 = sum(j2,1)';  
% [g,g2, abs(g-g2)]
% EstimOpt.NumGrad = 0;
% 
% pause;

if EstimOpt.NumGrad == 0
    if isequal(OptimOpt.Hessian,'user-supplied') == 1
        [LL, j] = LLfun(b0);
        f = -sum(Weights.*LL);
        j = j.*Weights(:, ones(size(j,2),1));
        g = -sum(j); ...
        h = j'*j;
    else
        [LL, j] = LLfun(b0);
        f = -sum(Weights.*LL);

        g = -sum(j.*Weights(:, ones(size(j,2),1))); ...
    end
else
    LL = LLfun(b0);
    f = -sum(Weights.*LL);
    if isequal(OptimOpt.Hessian,'user-supplied') == 1
        j = -numdiff(LLfun,LL,b0,isequal(OptimOpt.FinDiffType,'central'),ones(EstimOpt.NVarA,1));...
        j = j.*Weights(:, ones(size(j,2),1));    
        g = sum(j,1)'; ...
        h = j'*j;
    else
        j = -numdiff(LLfun,LL,b0,isequal(OptimOpt.FinDiffType,'central'),ones(EstimOpt.NVarA,1));...
        j = j.*Weights(:, ones(size(j,2),1));
        g = sum(j,1)';
    end 
end

