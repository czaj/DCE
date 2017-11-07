function [f,g,h] = LL_hlc_MATlike(YY,Xa,XXc,Xstr,Xmea,Xmea_exp,err_sliced,W,EstimOpt,OptimOpt,b0)

LLfun = @(B) LL_hlc(YY,Xa,XXc,Xstr,Xmea,Xmea_exp,err_sliced,EstimOpt,B);

LL = LLfun(b0);

if isequal(OptimOpt.GradObj,'on')
    if isequal(OptimOpt.Hessian,'user-supplied')
        j = numdiff(LLfun,LL,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
        j = j.*W;
        g = sum(j,1)';
        h = j'*j;
    else % OptimOpt.Hessian ~= 'user-supplied'
        j = numdiff(LLfun,LL,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
        j = j.*W;
        g = sum(j,1)';
    end
end
LL = LL.*W;
f = sum(LL);
