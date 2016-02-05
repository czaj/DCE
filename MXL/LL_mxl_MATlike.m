function [f,g,h]= LL_mxl_MATlike(YY,XXa,XXm,Xs,err,EstimOpt,OptimOpt,b0)

% save res_LL_mxl_MATlike;
% return

LLfun = @(B) LL_mxl(YY,XXa,XXm,Xs,err,EstimOpt,B);
[f,j] = LLfun(b0);
        j(:,EstimOpt.BActive ==0) = 0;
        g = sum(j,1)'; ...
EstimOpt.NumGrad = 1;
LLfun2 = @(B) LL_mxl(YY,XXa,XXm,Xs,err,EstimOpt,B);
f2 = LLfun2(b0);  
j2 = numdiff(LLfun2,f2,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);...
g2 = sum(j2,1)';  
[g,g2, abs(g-g2)]
EstimOpt.NumGrad = 0;
pause;

if isequal(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        if EstimOpt.ApproxHess == 1
            [f,j] = LLfun(b0);
            j(:,EstimOpt.BActive ==0) = 0;
            g = sum(j); ...
            if isequal(OptimOpt.Hessian,'user-supplied') == 1
                h = j'*j;
            end
        else
            [f,j,h] = LLfun(b0);
            j(:,EstimOpt.BActive ==0) = 0;
            g = sum(j); ...  
        end
    else % => EstimOpt.NumGrad == 1 
        f = LLfun(b0);
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            j = numdiff(LLfun,f,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
%             jacobian_new(@(B) LL_lcm(INPUT.YYcl,INPUT.Xa,INPUT.Xc,EstimOpt,B), b);
%             j = jacobian_new(LLfun, b0, 4, 1);...
            g = sum(j,1)'; ...
            h = j'*j;
        else % OptimOpt.Hessian ~= 'user-supplied'
%             j = numdiff_mxl(@(B0) LL_mxl(YY,XXa,XXm,XXs,TIMES,err_sliced,EstimOpt,B0),b0,EstimOpt);...
            j = numdiff(LLfun,f,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);...
            g = sum(j,1)';
        end
    end
else % No gradient
    EstimOpt.NumGrad = 1;
    f = LLfun(b0);   
end

f = sum(f);

end