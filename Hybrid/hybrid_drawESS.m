function [ess,ll] = hybrid_drawESS(model,args)
%HYBRID_DRAWESS Effective number of simulation draws per respondent.
%
% Syntax:   ess = hybrid_drawESS(model,args)
%           [ess,ll] = hybrid_drawESS(model,args)
%
% Inputs:
%    model - 'HMNL', 'HMXL', 'HLC' or 'MIMIC'
%    args  - the cell array the module passes to hybrid_likelihood, that is the
%            model's inputs followed by EstimOpt and the parameter vector
%
% Outputs:
%    ess   - NP x 1, Kish's effective sample size of the weights the draws enter
%            each respondent's simulated likelihood with:
%                L_i = (1/R) sum_r joint_ir,  w_ir = joint_ir/sum_r joint_ir,
%                ess_i = 1/sum_r w_ir^2.
%            In a hybrid model joint_ir is the product of the choice
%            probabilities and the likelihood of that respondent's indicators at
%            draw r, so the weights are the posterior of the latent variables
%            given everything the respondent told us. The more indicators, the
%            tighter that posterior is relative to the distribution the draws
%            come from, and the fewer draws carry the likelihood. ess_i = R means
%            every draw contributes equally; ess_i = 1 means one draw does.
%            For quasi-random draws it is conservative: it is the
%            independent-draw formula and ignores the balance of the net.
%    ll    - the simulated log-likelihood recomputed here, which the tests
%            compare with the module's own to show that this reads the same
%            integral the module maximises.
%
% It runs the same per-respondent blocks the likelihood does (hybrid_block),
% so there is no second copy of the model here, and it costs one likelihood
% evaluation without the gradient. Nothing is changed.

if exist('hybrid_likelihood','file') ~= 2
    error('DCE:Hybrid:Missing','the public hybrid path is not on the MATLAB path');
end
o = args{end-1};
b = args{end};
d = hybrid_prepare(model,args);
p = hybrid_parameters(d,b);

N = o.NP;
block = 32;
if isfield(o,'HybridBlockSize') && ~isempty(o.HybridBlockSize) && o.HybridBlockSize >= 1
    block = min(o.HybridBlockSize,N);
end
starts = 1:block:N;
ess = NaN(N,1);
f = NaN(N,1);
for j = 1:numel(starts)
    ids = starts(j):min(starts(j)+block-1,N);
    [fj,~,ej] = hybrid_block(d,p,ids,false);
    f(ids) = fj;
    ess(ids) = ej;
end
ll = -sum(f);
end
