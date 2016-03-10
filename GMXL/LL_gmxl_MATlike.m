function [f,g,h]= LL_gmxl_MATlike(YY,XXa,XXm,XXs,XXt,err_sliced,W,EstimOpt,OptimOpt,b0)

% save res_LL_mxl_MATlike;
% return

LLfun = @(B) LL_gmxl(YY,XXa,XXm,XXs,XXt,err_sliced,EstimOpt,B);

f = LLfun(b0);
% f = sum(LL_mxl(YY,XXa,XXm,XXs,TIMES,err_sliced,EstimOpt,b0));

if isequal(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        if EstimOpt.NVarS + EstimOpt.NVarM > 0
            error ('Error - analytical gradient not available for models with covariates of means or scale.')
        end
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            j = gmxl_gr8a(YY,XXa,TIMES,err_mtx,EstimOpt,b0); ... 
            j = j.*W(:, ones(1,size(j,2)));
            g = sum(j,1)'; ...
            h = j'*j;
        else
            g = sum(gmxl_gr8a(YY,XXa,TIMES,err_mtx,EstimOpt,b0));
        end
    else % => EstimOpt.NumGrad == 1         
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            j = numdiff(LLfun,f,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);...
            j = j.*W(:, ones(1,size(j,2)));
            g = sum(j,1)'; ...
            h = j'*j;
        else % OptimOpt.Hessian ~= 'user-supplied'
            j = numdiff(LLfun,f,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);...
            j = j.*W(:, ones(1,size(j,2)));
            g = sum(j,1)';
        end
    end
end
f = f.*W;
f = sum(f);
% save res_LL_mxl_MATlike;

end