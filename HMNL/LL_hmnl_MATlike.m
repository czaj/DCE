function [LL,g,h] = LL_hmnl_MATlike(Y,Xa,Xm,Xs,X_str,X_mea,Xmea_exp,err_sliced,W,EstimOpt,OptimOpt,b0)
% LL_hmnl_MATlike MATLAB-optimizer wrapper for LL_hmnl.
%
% Robustness additions:
%   - keeps a global record of the last/best finite evaluation;
%   - turns non-finite line-search probes into finite penalties when
%     EstimOpt.RealMin >= 2;
%   - reports whether the problem came from f, the Jacobian, or parameters;
%   - guards against absurd parameter values before LL_hmnl is called.

global HMNL_OPT_STATE

if nargin < 11
    error('LL_hmnl_MATlike:NotEnoughInputs','Not enough inputs.');
end

EstimOpt = local_defaults(EstimOpt);
b0 = b0(:);
W = W(:);

if isempty(HMNL_OPT_STATE) || ~isstruct(HMNL_OPT_STATE)
    HMNL_OPT_STATE = struct();
end
if ~isfield(HMNL_OPT_STATE,'eval_count') || isempty(HMNL_OPT_STATE.eval_count)
    HMNL_OPT_STATE.eval_count = 0;
end
HMNL_OPT_STATE.eval_count = HMNL_OPT_STATE.eval_count + 1;

npar = length(b0);
h = [];
g = zeros(npar,1);
LLfun = @(B) LL_hmnl(Y,Xa,Xm,Xs,X_str,X_mea,Xmea_exp,err_sliced,EstimOpt,B);
evalTimer = tic;
if EstimOpt.DebugOptimizationEval ~= 0 && (mod(HMNL_OPT_STATE.eval_count-1,EstimOpt.DebugEvalEvery) == 0)
    local_progress_msg(sprintf('HMNL eval %d started: max|B|=%.4g, NRep=%d.',HMNL_OPT_STATE.eval_count,max(abs(b0)),EstimOpt.NRep));
end

% Hard parameter sanity check before expensive likelihood evaluation.
[maxAbsB,idxMaxB] = max(abs(b0));
if any(~isfinite(b0)) || maxAbsB > EstimOpt.MaxAbsB
    [LL,g,h] = local_parameter_penalty(b0,EstimOpt,OptimOpt);
    local_register_eval(b0,LL,g,[],[],EstimOpt,'parameter_limit');
    if EstimOpt.ReportOptimizationDiagnostics ~= 0
        local_warn_once(sprintf('HMNL objective penalty: parameter %d reached %.4g (limit %.4g).',idxMaxB,maxAbsB,EstimOpt.MaxAbsB),EstimOpt);
    end
    return
end

try
    if isequal(OptimOpt.GradObj,'on')
        if EstimOpt.NumGrad == 0
            [f,j] = LLfun(b0);
        else
            f = LLfun(b0);
            j = numdiff(LLfun,f,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);
        end

        if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive)
            EstimOpt.BActive = ones(1,npar);
        end
        EstimOpt.BActive = EstimOpt.BActive(:)';
        if length(EstimOpt.BActive) ~= npar
            error('LL_hmnl_MATlike:BActiveSize','EstimOpt.BActive has incorrect length.');
        end

        j(:,EstimOpt.BActive == 0) = 0;
        [f,j,diagInfo] = local_repair_fj(f,j,b0,EstimOpt);

        jw = j.*W;
        g = sum(jw,1)';
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            h = jw'*jw;
        end
    else
        EstimOpt.NumGrad = 1;
        f = LLfun(b0);
        j = [];
        [f,~,diagInfo] = local_repair_fj(f,[],b0,EstimOpt);
    end
catch ME
    if EstimOpt.RealMin >= 2
        [LL,g,h] = local_failure_penalty(b0,EstimOpt,OptimOpt,ME);
        local_register_eval(b0,LL,g,[],[],EstimOpt,'caught_error');
        if EstimOpt.ReportOptimizationDiagnostics ~= 0
            local_warn_once(['HMNL objective penalty after caught error: ' ME.message],EstimOpt);
        end
        return
    else
        rethrow(ME)
    end
end

f = f(:).*W;
LL = sum(f);

if ~isfinite(LL) || any(~isfinite(g))
    if EstimOpt.RealMin >= 2
        if ~isfinite(LL)
            bad = ~isfinite(f);
            f(bad) = EstimOpt.BadEvalPenalty/length(f);
            LL = sum(f);
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
if EstimOpt.DebugOptimizationEval ~= 0 && (mod(HMNL_OPT_STATE.eval_count-1,EstimOpt.DebugEvalEvery) == 0)
    local_progress_msg(sprintf('HMNL eval %d finished in %.2fs: LL=%.6g, max|g|=%.4g.',HMNL_OPT_STATE.eval_count,diagInfo.elapsed_seconds,LL,max(abs(g))));
end
local_register_eval(b0,LL,g,f,diagInfo,EstimOpt,'normal');

end

function EstimOpt = local_defaults(EstimOpt)
if ~isfield(EstimOpt,'RealMin') || isempty(EstimOpt.RealMin), EstimOpt.RealMin = 0; end
if ~isfield(EstimOpt,'ProbMin') || isempty(EstimOpt.ProbMin), EstimOpt.ProbMin = 1e-300; end
EstimOpt.ProbMin = max(realmin,min(EstimOpt.ProbMin,1e-6));
if ~isfield(EstimOpt,'MaxExp') || isempty(EstimOpt.MaxExp), EstimOpt.MaxExp = 50; end
if ~isfield(EstimOpt,'SigmaMin') || isempty(EstimOpt.SigmaMin), EstimOpt.SigmaMin = 1e-8; end
if ~isfield(EstimOpt,'MaxAbsB') || isempty(EstimOpt.MaxAbsB), EstimOpt.MaxAbsB = 1000; end
if ~isfield(EstimOpt,'MaxAbsGrad') || isempty(EstimOpt.MaxAbsGrad), EstimOpt.MaxAbsGrad = 1e12; end
if ~isfield(EstimOpt,'BadEvalPenalty') || isempty(EstimOpt.BadEvalPenalty), EstimOpt.BadEvalPenalty = 1e50; end
if ~isfield(EstimOpt,'PenaltySlope') || isempty(EstimOpt.PenaltySlope), EstimOpt.PenaltySlope = 1e6; end
if ~isfield(EstimOpt,'ReportOptimizationDiagnostics') || isempty(EstimOpt.ReportOptimizationDiagnostics), EstimOpt.ReportOptimizationDiagnostics = 1; end
if ~isfield(EstimOpt,'MaxDiagnosticWarnings') || isempty(EstimOpt.MaxDiagnosticWarnings), EstimOpt.MaxDiagnosticWarnings = 10; end
if ~isfield(EstimOpt,'DebugOptimizationEval') || isempty(EstimOpt.DebugOptimizationEval), EstimOpt.DebugOptimizationEval = 0; end
if ~isfield(EstimOpt,'DebugEvalEvery') || isempty(EstimOpt.DebugEvalEvery), EstimOpt.DebugEvalEvery = 1; end
end

function [f,j,diagInfo] = local_repair_fj(f,j,b,EstimOpt)
f = f(:);
diagInfo = struct();
diagInfo.eval_time = datestr(now,'yyyy-mm-dd HH:MM:SS');
diagInfo.max_abs_B = max(abs(b));
diagInfo.n_nonfinite_f = sum(~isfinite(f));
diagInfo.max_f = max(f(isfinite(f)));
diagInfo.min_f = min(f(isfinite(f)));
if isempty(diagInfo.max_f), diagInfo.max_f = NaN; end
if isempty(diagInfo.min_f), diagInfo.min_f = NaN; end

if ~isempty(j)
    badJ = ~isfinite(j);
    diagInfo.n_nonfinite_j = sum(badJ(:));
    if EstimOpt.ReportOptimizationDiagnostics ~= 0 || EstimOpt.RealMin >= 3
        jtmp = j;
        jtmp(badJ) = 0;
        diagInfo.max_abs_j = max(abs(jtmp(:)));
    else
        diagInfo.max_abs_j = NaN;
    end
else
    badJ = [];
    diagInfo.n_nonfinite_j = 0;
    diagInfo.max_abs_j = NaN;
end

if EstimOpt.RealMin >= 2
    if diagInfo.n_nonfinite_f > 0
        f(~isfinite(f)) = EstimOpt.BadEvalPenalty/length(f);
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

function [LL,g,h] = local_parameter_penalty(b,EstimOpt,OptimOpt)
npar = length(b);
excess = max(0,abs(b) - EstimOpt.MaxAbsB);
LL = EstimOpt.BadEvalPenalty + EstimOpt.PenaltySlope*sum(excess.^2);
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

function [LL,g,h] = local_failure_penalty(b,EstimOpt,OptimOpt,ME)
[LL,g,h] = local_parameter_penalty(b,EstimOpt,OptimOpt);
if isfield(ME,'message') && ~isempty(ME.message)
    LL = LL + min(1e20,length(ME.message));
end
end

function local_register_eval(b,LL,g,f,diagInfo,EstimOpt,status)
global HMNL_OPT_STATE
HMNL_OPT_STATE.last_B = b(:);
HMNL_OPT_STATE.last_LL = LL;
HMNL_OPT_STATE.last_g = g;
HMNL_OPT_STATE.last_status = status;
HMNL_OPT_STATE.last_diag = diagInfo;
if isfinite(LL) && all(isfinite(g))
    if ~isfield(HMNL_OPT_STATE,'best_LL') || isempty(HMNL_OPT_STATE.best_LL) || LL < HMNL_OPT_STATE.best_LL
        HMNL_OPT_STATE.best_LL = LL;
        HMNL_OPT_STATE.best_B = b(:);
        HMNL_OPT_STATE.best_g = g;
        HMNL_OPT_STATE.best_status = status;
        HMNL_OPT_STATE.best_diag = diagInfo;
    end
end

if EstimOpt.ReportOptimizationDiagnostics ~= 0
    bad = false;
    if isstruct(diagInfo)
        bad = (isfield(diagInfo,'n_nonfinite_f') && diagInfo.n_nonfinite_f > 0) || ...
              (isfield(diagInfo,'n_nonfinite_j') && diagInfo.n_nonfinite_j > 0) || ...
              (isfield(diagInfo,'max_abs_j') && isfinite(diagInfo.max_abs_j) && diagInfo.max_abs_j > EstimOpt.MaxAbsGrad);
    end
    if bad
        msg = sprintf('HMNL diagnostic: status=%s, LL=%.6g, nonfinite f=%d, nonfinite jac=%d, max|jac|=%.4g, max|B|=%.4g.', ...
            status,LL,diagInfo.n_nonfinite_f,diagInfo.n_nonfinite_j,diagInfo.max_abs_j,diagInfo.max_abs_B);
        local_warn_once(msg,EstimOpt);
    end
end
end

function local_progress_msg(msg)
try
    disp(msg);
catch
end
end

function local_warn_once(msg,EstimOpt)
global HMNL_OPT_STATE
if ~isfield(HMNL_OPT_STATE,'warning_count') || isempty(HMNL_OPT_STATE.warning_count)
    HMNL_OPT_STATE.warning_count = 0;
end
if HMNL_OPT_STATE.warning_count < EstimOpt.MaxDiagnosticWarnings
    HMNL_OPT_STATE.warning_count = HMNL_OPT_STATE.warning_count + 1;
    try
        fprintf(2,'WARNING: %s\n',msg);
    catch
        warning('HMNL:OptimizationDiagnostic','%s',msg);
    end
end
end
