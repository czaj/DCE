function [f,g,h]= LL_mxl_MATlike(YY,XXa,XXm,Xs,err,W,EstimOpt,OptimOpt,b0)

% save res_LL_mxl_MATlike;
% return

LLfun = @(B) LL_mxl(YY,XXa,XXm,Xs,err,EstimOpt,B);

if isequal(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        if EstimOpt.ApproxHess == 1
            [f,j] = LLfun(b0);
            j(:,EstimOpt.BActive ==0) = 0;
            j = j.*W(:, ones(1,size(j,2)));
            g = sum(j);
            if isequal(OptimOpt.Hessian,'user-supplied') == 1
                h = j'*j;
            end
        else
            [f,j,h] = LLfun(b0);
            j(:,EstimOpt.BActive ==0) = 0;
            j = j.*W(:, ones(1,size(j,2)));
            g = sum(j);
        end
    else % => EstimOpt.NumGrad == 1
        f = LLfun(b0);
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            j = numdiff(LLfun,f,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
            %             jacobian_new(@(B) LL_lcm(INPUT.YYcl,INPUT.Xa,INPUT.Xc,EstimOpt,B), b);
            %             j = jacobian_new(LLfun, b0, 4, 1);...
            j = j.*W(:, ones(1,size(j,2)));
            g = sum(j,1)';
            h = j'*j;
        else % OptimOpt.Hessian ~= 'user-supplied'
            %             j = numdiff_mxl(@(B0) LL_mxl(YY,XXa,XXm,XXs,TIMES,err_sliced,EstimOpt,B0),b0,EstimOpt);...
            j = numdiff(LLfun,f,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);...
                j = j.*W(:, ones(1,size(j,2)));
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