function [INPUT, EstimOpt, OptimOpt] = setDisplayOptions(INPUT, EstimOpt, OptimOpt)
% Needs to be refactored
% Possibly couple of instructions should be moved to setOptimizationOptions

if ((isfield(EstimOpt,'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || ...
        sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    OptimOpt.GradObj = 'on';
end

% if EstimOpt.NVarS > 0 && EstimOpt.NumGrad == 0 && any(isnan(INPUT.Xa(:)))
% 	EstimOpt.NumGrad = 1;
%     OptimOpt.GradObj = 'off';
% 	cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient off - covariates of scale not supported by analytical gradient \n')
% end

% TODO
% jezeli zostanie zrobiony gradient analityczny to do zmiany
if (any(EstimOpt.Dist > 1) && ~any(EstimOpt.Dist == 4)) && EstimOpt.NumGrad == 0
    EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'off';
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient off - analytical gradient available for normally, lognormally or Weibull distributed parameters only \n')
end


% This is not necessary any more (?)
% if any(var(EstimOpt.NAltMissInd)) ~= 0 && EstimOpt.NumGrad == 0 && EstimOpt.WTP_space > 1
%     EstimOpt.NumGrad = 1;
%     OptimOpt.GradObj = 'off';
%     cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient off - analytical gradient not available when the number of alternatives differs and WTP_space > 1  \n')
% end

if ((isfield(EstimOpt,'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || ...
        sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'on';
end

% if EstimOpt.NVarNLT > 0 && EstimOpt.NLTType == 2 && EstimOpt.NumGrad == 0
% 	EstimOpt.NumGrad = 1;
% 	if EstimOpt.Display ~= 0
%         cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient to numerical - Yeo-Johnston transformation not supported by analytical gradient \n')
% 	end
% end

if (isfield(EstimOpt,'ConstVarActive') == 0 || EstimOpt.ConstVarActive == 0) && ...
        isequal(OptimOpt.Algorithm,'quasi-newton') && isequal(OptimOpt.Hessian,'user-supplied')
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied Hessian off - quasi-newton algorithm does not use it anyway \n')
    OptimOpt.Hessian = 'off';
end

if EstimOpt.NumGrad == 1 && EstimOpt.ApproxHess == 0
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied exact Hessian off - exact Hessian only available if analythical gradient on \n')
    EstimOpt.ApproxHess = 1;
end

if EstimOpt.NVarS > 0 && (EstimOpt.ApproxHess == 0 || EstimOpt.HessEstFix == 4)
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied exact Hessian off - exact Hessian not available for models with covariates of scale \n')
    EstimOpt.ApproxHess = 1;
    EstimOpt.HessEstFix = 0;
end

if EstimOpt.NVarM > 0 && (EstimOpt.ApproxHess == 0 || EstimOpt.HessEstFix == 4)
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied exact Hessian off - exact Hessian not available for models with covariates of means \n')
    EstimOpt.ApproxHess = 1;
    EstimOpt.HessEstFix = 0;
end

if EstimOpt.WTP_space > 0 && (EstimOpt.ApproxHess == 0 || EstimOpt.HessEstFix == 4)
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied exact Hessian off - exact Hessian not available for models in WTP-space \n')
    EstimOpt.ApproxHess = 1;
    EstimOpt.HessEstFix = 0;
end

if any(isnan(INPUT.Xa(:))) == 1 && (EstimOpt.ApproxHess == 0 || EstimOpt.HessEstFix == 4)
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied exact Hessian off - exact Hessian not available with missing data \n')
    EstimOpt.ApproxHess = 1;
    EstimOpt.HessEstFix = 0;
end

if any(EstimOpt.Dist ~= 0) && (EstimOpt.ApproxHess == 0 || EstimOpt.HessEstFix == 4)
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied exact Hessian off - exact Hessian available for models with normally distributed parameters only \n')
    EstimOpt.ApproxHess = 1;
    EstimOpt.HessEstFix = 0;
end

if EstimOpt.NVarNLT > 0 && (EstimOpt.ApproxHess == 0 || EstimOpt.HessEstFix == 4)
    cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied exact Hessian off - exact Hessian not available for models with non-linear transformation(s) of variable(s) \n')
    EstimOpt.ApproxHess = 1;
    EstimOpt.HessEstFix = 0;
end

if any(INPUT.W ~= 1) && ...
        ((EstimOpt.ApproxHess == 0 && EstimOpt.NumGrad == 0) || EstimOpt.HessEstFix == 4)
    INPUT.W = ones(EstimOpt.NP,1);
    cprintf(rgb('DarkOrange'),'WARNING: Setting all weights to 1, they are not supported with analytical hessian \n')
end

if EstimOpt.RobustStd == 1 && (EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2)
    EstimOpt.RobustStd = 0;
    cprintf(rgb('DarkOrange'),'WARNING: Setting off robust standard errors, they do not matter for BHHH aproximation of hessian \n')
end

if  any(EstimOpt.Dist >= 3 & EstimOpt.Dist <= 5) && EstimOpt.NVarM ~= 0
    error('Covariates of means do not work with triangular/weibull/sinh-arcsinh distributions')
end

fprintf('\n')
cprintf('Optimization algorithm: '); cprintf('*Black',[OptimOpt.Algorithm '\n'])

if strcmp(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        cprintf('Gradient: '); cprintf('*Black','user-supplied, analytical \n')
    else
        cprintf('Gradient: '); cprintf('*Black',['user-supplied, numerical, ' OptimOpt.FinDiffType '\n'])
    end
else
    cprintf('Gradient: '); cprintf('*Black',['built-in, ' OptimOpt.FinDiffType '\n'])
end

if isequal(OptimOpt.Algorithm,'quasi-newton')
    cprintf('Hessian: '); cprintf('*Black','off, ')
    switch EstimOpt.HessEstFix
        case 0
            cprintf('*Black','retained from optimization \n')
        case 1
            cprintf('*Black','ex-post calculated using BHHH \n')
        case 2
            cprintf('*Black','ex-post calculated using high-precision BHHH \n')
        case 3
            cprintf('*Black','ex-post calculated numerically \n')
        case 4
            cprintf('*Black','ex-post calculated analytically \n')
    end
else
    if strcmp(OptimOpt.Hessian,'user-supplied')
        if EstimOpt.ApproxHess == 1
            cprintf('Hessian: '); cprintf('*Black','user-supplied, BHHH, ')
        else
            cprintf('Hessian: '); cprintf('*Black','user-supplied, analytical, ')
        end
    else
        cprintf('Hessian: '); cprintf('*Black',['built-in, ' OptimOpt.HessUpdate ', '])
    end
    switch EstimOpt.HessEstFix
        case 0
            cprintf('*Black','retained from optimization \n')
        case 1
            cprintf('*Black','ex-post calculated using BHHH \n')
        case 2
            cprintf('*Black','ex-post calculated using high-precision BHHH \n')
        case 3
            cprintf('*Black','ex-post calculated numerically \n')
        case 4
            cprintf('*Black','ex-post calculated analytically \n')
    end
end
fprintf('\n')

end

