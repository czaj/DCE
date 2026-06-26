function [f,g,h]= LL_mxl_MATlike(YY,XXa,XXm,Xs,err,W,EstimOpt,OptimOpt,b0)
% Function calculating loglikelihood function value (f), its gradient (g)
% and hessian (h).

global MXL_OPT_STATE

EstimOpt = local_defaults(EstimOpt);
b0 = b0(:);
W = W(:);
npar = length(b0);
h = [];
g = zeros(npar,1);

if isempty(MXL_OPT_STATE) || ~isstruct(MXL_OPT_STATE)
    MXL_OPT_STATE = struct();
end
if ~isfield(MXL_OPT_STATE,'eval_count') || isempty(MXL_OPT_STATE.eval_count)
    MXL_OPT_STATE.eval_count = 0;
end
MXL_OPT_STATE.eval_count = MXL_OPT_STATE.eval_count + 1;

LLfun = @(B) LL_mxl(YY,XXa,XXm,Xs,err,EstimOpt,B);
evalTimer = tic;

[maxAbsB,idxMaxB] = max(abs(b0));
if any(~isfinite(b0)) || maxAbsB > EstimOpt.MaxAbsB
    [f,g,h] = local_parameter_penalty(b0,EstimOpt,OptimOpt);
    local_register_eval(b0,f,g,[],[],EstimOpt,'parameter_limit');
    if EstimOpt.ReportOptimizationDiagnostics ~= 0
        local_warn_once(sprintf('MXL objective penalty: parameter %d reached %.4g (limit %.4g).',idxMaxB,maxAbsB,EstimOpt.MaxAbsB),EstimOpt);
    end
    return
end

try
    if isequal(OptimOpt.GradObj,'on')
        if EstimOpt.NumGrad == 0
            if EstimOpt.ApproxHess == 1
                [fv,j] = LLfun(b0);
            else
                [fv,j,hraw] = LLfun(b0);
                h = hraw;
            end
        else
            fv = LLfun(b0);
            j = numdiff(LLfun,fv,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
        end

        if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive)
            EstimOpt.BActive = ones(1,npar);
        end
        EstimOpt.BActive = EstimOpt.BActive(:)';
        if length(EstimOpt.BActive) ~= npar
            error('LL_mxl_MATlike:BActiveSize','EstimOpt.BActive has incorrect length.');
        end

        j(:,EstimOpt.BActive == 0) = 0;
        [fv,j,diagInfo] = local_repair_fj(fv,j,b0,EstimOpt);
        jw = j.*W;
        g = sum(jw,1)';
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            if EstimOpt.ApproxHess == 1 || isempty(h) || any(~isfinite(h(:)))
                h = jw'*jw;
            end
        end
    else
        EstimOpt.NumGrad = 1;
        fv = LLfun(b0);
        [fv,~,diagInfo] = local_repair_fj(fv,[],b0,EstimOpt);
    end
catch ME
    if EstimOpt.RealMin >= 2
        [f,g,h] = local_failure_penalty(b0,EstimOpt,OptimOpt,ME);
        local_register_eval(b0,f,g,[],[],EstimOpt,'caught_error');
        if EstimOpt.ReportOptimizationDiagnostics ~= 0
            local_warn_once(['MXL objective penalty after caught error: ' ME.message],EstimOpt);
        end
        return
    else
        rethrow(ME)
    end
end

fv = fv(:).*W;
f = sum(fv); % -simulated log likelihood

if ~isfinite(f) || any(~isfinite(g))
    if EstimOpt.RealMin >= 2
        if ~isfinite(f)
            bad = ~isfinite(fv);
            fv(bad) = EstimOpt.BadEvalPenalty/length(fv);
            f = sum(fv);
        end
        if any(~isfinite(g))
            g(~isfinite(g)) = 0;
        end
        if EstimOpt.RealMin >= 3 && any(abs(g) > EstimOpt.MaxAbsGrad)
            g = max(-EstimOpt.MaxAbsGrad,min(EstimOpt.MaxAbsGrad,g));
        end
        if isequal(OptimOpt.Hessian,'user-supplied') == 1 && (isempty(h) || any(~isfinite(h(:))))
            h = eye(npar);
        end
    end
end

diagInfo.elapsed_seconds = toc(evalTimer);
local_register_eval(b0,f,g,fv,diagInfo,EstimOpt,'normal');
end

function EstimOpt = local_defaults(EstimOpt)
if ~isfield(EstimOpt,'RealMin') || isempty(EstimOpt.RealMin), EstimOpt.RealMin = 0; end
if ~isfield(EstimOpt,'NumGrad') || isempty(EstimOpt.NumGrad), EstimOpt.NumGrad = 0; end
if ~isfield(EstimOpt,'ApproxHess') || isempty(EstimOpt.ApproxHess), EstimOpt.ApproxHess = 1; end
if ~isfield(EstimOpt,'MaxAbsB') || isempty(EstimOpt.MaxAbsB), EstimOpt.MaxAbsB = 1000; end
if ~isfield(EstimOpt,'MaxAbsGrad') || isempty(EstimOpt.MaxAbsGrad), EstimOpt.MaxAbsGrad = 1e12; end
if ~isfield(EstimOpt,'BadEvalPenalty') || isempty(EstimOpt.BadEvalPenalty), EstimOpt.BadEvalPenalty = 1e50; end
if ~isfield(EstimOpt,'PenaltySlope') || isempty(EstimOpt.PenaltySlope), EstimOpt.PenaltySlope = 1e6; end
if ~isfield(EstimOpt,'ReportOptimizationDiagnostics') || isempty(EstimOpt.ReportOptimizationDiagnostics), EstimOpt.ReportOptimizationDiagnostics = 1; end
if ~isfield(EstimOpt,'MaxDiagnosticWarnings') || isempty(EstimOpt.MaxDiagnosticWarnings), EstimOpt.MaxDiagnosticWarnings = 10; end
end

function [fv,j,diagInfo] = local_repair_fj(fv,j,b,EstimOpt)
fv = fv(:);
diagInfo = struct();
diagInfo.eval_time = datestr(now,'yyyy-mm-dd HH:MM:SS');
diagInfo.max_abs_B = max(abs(b));
diagInfo.n_nonfinite_f = sum(~isfinite(fv));
diagInfo.max_f = max(fv(isfinite(fv)));
diagInfo.min_f = min(fv(isfinite(fv)));
if isempty(diagInfo.max_f), diagInfo.max_f = NaN; end
if isempty(diagInfo.min_f), diagInfo.min_f = NaN; end

if ~isempty(j)
    badJ = ~isfinite(j);
    diagInfo.n_nonfinite_j = sum(badJ(:));
    jtmp = j;
    jtmp(badJ) = 0;
    diagInfo.max_abs_j = max(abs(jtmp(:)));
else
    badJ = [];
    diagInfo.n_nonfinite_j = 0;
    diagInfo.max_abs_j = NaN;
end

if EstimOpt.RealMin >= 2
    if diagInfo.n_nonfinite_f > 0
        fv(~isfinite(fv)) = EstimOpt.BadEvalPenalty/length(fv);
    end
    if ~isempty(j)
        if diagInfo.n_nonfinite_j > 0
            j(badJ) = 0;
        end
        if EstimOpt.RealMin >= 3
            j = max(-EstimOpt.MaxAbsGrad,min(EstimOpt.MaxAbsGrad,j));
        end
    end
end
end

function [f,g,h] = local_parameter_penalty(b,EstimOpt,OptimOpt)
npar = length(b);
excess = max(0,abs(b) - EstimOpt.MaxAbsB);
f = EstimOpt.BadEvalPenalty + EstimOpt.PenaltySlope*sum(excess.^2);
g = 2*EstimOpt.PenaltySlope*excess.*sign(b);
g(~isfinite(g)) = 0;
if ~any(excess)
    g = zeros(npar,1);
end
if isequal(OptimOpt.Hessian,'user-supplied') == 1
    h = eye(npar)*EstimOpt.PenaltySlope;
else
    h = [];
end
end

function [f,g,h] = local_failure_penalty(b,EstimOpt,OptimOpt,ME)
[f,g,h] = local_parameter_penalty(b,EstimOpt,OptimOpt);
if isfield(ME,'message') && ~isempty(ME.message)
    f = f + min(1e20,length(ME.message));
end
end

function local_register_eval(b,f,g,fv,diagInfo,EstimOpt,status)
global MXL_OPT_STATE
MXL_OPT_STATE.last_B = b(:);
MXL_OPT_STATE.last_LL = f;
MXL_OPT_STATE.last_g = g;
MXL_OPT_STATE.last_status = status;
MXL_OPT_STATE.last_f = fv;
MXL_OPT_STATE.last_diag = diagInfo;
if isfinite(f) && all(isfinite(g))
    if ~isfield(MXL_OPT_STATE,'best_LL') || isempty(MXL_OPT_STATE.best_LL) || f < MXL_OPT_STATE.best_LL
        MXL_OPT_STATE.best_LL = f;
        MXL_OPT_STATE.best_B = b(:);
        MXL_OPT_STATE.best_g = g;
        MXL_OPT_STATE.best_status = status;
        MXL_OPT_STATE.best_diag = diagInfo;
    end
end

if EstimOpt.ReportOptimizationDiagnostics ~= 0 && isstruct(diagInfo)
    bad = (isfield(diagInfo,'n_nonfinite_f') && diagInfo.n_nonfinite_f > 0) || ...
          (isfield(diagInfo,'n_nonfinite_j') && diagInfo.n_nonfinite_j > 0) || ...
          (isfield(diagInfo,'max_abs_j') && isfinite(diagInfo.max_abs_j) && diagInfo.max_abs_j > EstimOpt.MaxAbsGrad);
    if bad
        msg = sprintf('MXL diagnostic: status=%s, LL=%.6g, nonfinite f=%d, nonfinite jac=%d, max|jac|=%.4g, max|B|=%.4g.', ...
            status,f,diagInfo.n_nonfinite_f,diagInfo.n_nonfinite_j,diagInfo.max_abs_j,diagInfo.max_abs_B);
        local_warn_once(msg,EstimOpt);
    end
end
end

function local_warn_once(msg,EstimOpt)
global MXL_OPT_STATE
if ~isfield(MXL_OPT_STATE,'warning_count') || isempty(MXL_OPT_STATE.warning_count)
    MXL_OPT_STATE.warning_count = 0;
end
if MXL_OPT_STATE.warning_count < EstimOpt.MaxDiagnosticWarnings
    MXL_OPT_STATE.warning_count = MXL_OPT_STATE.warning_count + 1;
    try
        fprintf(2,'WARNING: %s\n',msg);
    catch
        warning('MXL:OptimizationDiagnostic','%s',msg);
    end
end
end
