function [ess,ll] = mxl_drawESS(INPUT,EstimOpt,err,b0)
%MXL_DRAWESS Effective number of simulation draws per respondent in MXL.
%
% Syntax:   ess = mxl_drawESS(INPUT,EstimOpt,err,b0)
%           [ess,ll] = mxl_drawESS(INPUT,EstimOpt,err,b0)
%
% Inputs:
%    INPUT   - the module's INPUT after reshaping (INPUT.XXa and INPUT.YY)
%    EstimOpt - the module's EstimOpt
%    err     - the module's draws, NVarA x NRep*NP with NRep consecutive
%              columns per respondent
%    b0      - the parameter vector the diagnostics are computed at
%
% Outputs:
%    ess     - NP x 1, Kish's effective sample size of the weights the draws
%              enter each respondent's simulated likelihood with:
%                  L_i = (1/R) sum_r P_i(beta_ir),  w_ir = P_i(beta_ir)/sum_r,
%                  ess_i = 1/sum_r w_ir^2 = (sum_r P_ir)^2/sum_r P_ir^2.
%              ess_i = R means every draw contributes equally, ess_i = 1 means
%              a single draw carries all the weight. This statistic measures
%              weight concentration for both random and quasi-random draws.
%    ll      - the simulated log-likelihood recomputed here. Used by the tests
%              to prove that this function integrates the same thing LL_mxl
%              does; it is not meant to replace the module's own value.
%
% Supported specifications: normal and lognormal random parameters, diagonal
% or full covariance, preference space or WTP space, whole choice tasks
% missing. Everything else (other distributions, covariates of means or scale,
% transformed variables, alternative-level missing data) raises an error, which
% the caller reports as the diagnostic being unavailable. The loop runs one
% respondent at a time, so it needs no NP x NRep workspace.

NAlt = EstimOpt.NAlt;
NCT = EstimOpt.NCT;
NP = EstimOpt.NP;
NRep = EstimOpt.NRep;
NVarA = EstimOpt.NVarA;
Dist = EstimOpt.Dist(:)';
WTP_space = EstimOpt.WTP_space;

if ~all(Dist == 0 | Dist == 1 | Dist == -1)
    error('mxl_drawESS:Distribution', ...
        'only normal and lognormal random parameters are supported (EstimOpt.Dist)');
end
if optField(EstimOpt,'NVarM',0) ~= 0 || optField(EstimOpt,'NVarS',0) ~= 0 || ...
        optField(EstimOpt,'NVarNLT',0) ~= 0 || optField(EstimOpt,'Johnson',0) ~= 0
    error('mxl_drawESS:Specification', ...
        'covariates of means or scale, transformed variables and Johnson distributions are not supported');
end
if WTP_space > 1
    error('mxl_drawESS:WTPspace','only one cost coefficient is supported in WTP space');
end
if size(err,2) ~= NRep*NP || size(err,1) ~= NVarA
    error('mxl_drawESS:Draws','the draw matrix is not NVarA x NRep*NP');
end

b0 = b0(:);
ExpB = optField(EstimOpt,'ExpB',[]);
if ~isempty(ExpB)
    b0(ExpB) = exp(b0(ExpB));
end
b0a = b0(1:NVarA);
if EstimOpt.FullCov == 0
    VC = diag(b0(NVarA+1:NVarA*2));
else
    VC = tril(ones(NVarA));
    VC(VC == 1) = b0(NVarA+1:NVarA+sum(1:NVarA));
end

XXa = INPUT.XXa;
YY = INPUT.YY;
idx = reshape(1:NAlt*NCT,[NAlt,NCT]);
ess = NaN(NP,1);
lli = zeros(NP,1); % a panel without observed choices contributes log(1)

for n = 1:NP
    Xn = XXa(:,:,n);
    yn = YY(:,n);
    ok = ~isnan(yn) & all(~isnan(Xn),2);
    okm = reshape(ok,[NAlt,NCT]);
    if isfield(EstimOpt,'MissingCT') && ~isempty(EstimOpt.MissingCT)
        okm(:,logical(EstimOpt.MissingCT(:,n))) = false;
    end
    keep = all(okm,1);
    if any(any(okm(:,~keep)))
        error('mxl_drawESS:MissingData', ...
            'respondent %d has a choice task with only some alternatives present',n);
    end
    nt = sum(keep);
    if nt == 0
        continue
    end
    rows = reshape(idx(:,keep),[],1);

    B = b0a + VC*err(:,(n-1)*NRep+1:n*NRep);     % NVarA x NRep
    B(Dist == 1,:) = exp(B(Dist == 1,:));
    if WTP_space == 1
        costs = optField(EstimOpt,'WTP_matrix',NVarA);
        B(1:NVarA-1,:) = B(1:NVarA-1,:).*B(costs,:);
    end

    logP = mxl_log_panel(Xn(rows,:),yn(rows)==1,B,NAlt);
    m = max(logP);
    if any(isnan(logP)) || ~isfinite(m)
        error('mxl_drawESS:Nonfinite','No finite log-panel integral for respondent %d.',n);
    end
    w = exp(logP-m);
    ess(n) = sum(w)^2/sum(w.^2);
    lli(n) = m + log(sum(w)/NRep);
end

ll = sum(lli);
end


function v = optField(S,name,default)
if isstruct(S) && isfield(S,name) && ~isempty(S.(name))
    v = S.(name);
else
    v = default;
end
end
