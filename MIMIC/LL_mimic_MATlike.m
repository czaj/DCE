function [LL,g,h] = LL_mimic_MATlike(Xstr,Xmea,Xmea_exp,err,W,EstimOpt,OptimOpt,b0)

LLfun = @(B) LL_mimic(Xstr,Xmea,Xmea_exp,err,EstimOpt,B);

if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive)
    EstimOpt.BActive = ones(1,length(b0));
else
    EstimOpt.BActive = EstimOpt.BActive(:)';
end

if isequal(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        [f,j] = LLfun(b0);
        j(:,EstimOpt.BActive == 0) = 0;
        j = j.*W;
        g = sum(j,1)';
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            h = j'*j;
        end
    else % => EstimOpt.NUMGRAD == 1 
        f = LLfun(b0);  
        j = numdiff(LLfun,f,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
        j = j.*W;
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

