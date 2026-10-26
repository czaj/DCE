function test_hybrid_optimization
% Coordinate chain rules, fixed nonzero values, callbacks, Hessian and cutoff delta.
Q = [4 1 0;1 3 0;0 0 2]; target = [2;-1;7]; initial = [-4;5;7]; active = [1;1;0];
o = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,...
    'Display','off','OptimalityTolerance',1e-9,'StepTolerance',1e-12,'OutputFcn',@checkCallback);
seen = false;
for whiten = 0:1
    H0 = []; if whiten, H0 = Q; end
    [b,f,flag,out,g,H] = hybrid_optimize(@quadratic,initial,o,active,H0);
    assert(flag > 0 && norm(b-target) < 1e-6 && f < 1e-12 && b(3) == 7);
    assert(norm(g) < 1e-6 && norm(H(1:2,1:2)-Q(1:2,1:2),'fro') < 1e-5);
    assert(all(H(3,:) == 0) && out.hybridActiveParameters == 2 && seen);
end
tr = o; tr.Algorithm = 'trust-region'; tr.HessianFcn = 'objective';
[b,~,flag] = hybrid_optimize(@quadratic,initial,tr,active,Q);
assert(flag > 0 && norm(b-target) < 1e-6);
num = o; num.SpecifyObjectiveGradient = false;
num.TypicalX = [1;2;3]; num.FiniteDifferenceStepSize = [1e-7;2e-7;3e-7];
[b,~,flag] = hybrid_optimize(@quadratic,initial,num,active);
assert(flag > 0 && norm(b-target) < 1e-4 && b(3) == 7);
stop = o; stop.OutputFcn = {@checkCallback,@(~,~,state) strcmp(state,'iter')};
[~,~,flag] = hybrid_optimize(@quadratic,initial,stop,active);
assert(flag == -1,'User stop requests must survive the coordinate map.');
bad = Q; bad(1,1) = -4; caught = false;
try, hybrid_optimize(@quadratic,initial,o,active,bad); catch ex, caught = strcmp(ex.identifier,'DCE:Hybrid:StartHessian'); end
assert(caught,'Indefinite curvature must not be silently repaired.');
raw = [-2;log(.7);log(1.3)]; V = [.4 .02 -.01;.02 .2 .03;-.01 .03 .3];
[threshold,se,p] = hybrid_cutoffs(raw,V);
J = [1 0 0;1 .7 0;1 .7 1.3];
assert(norm(threshold-[-2;-1.3;0]) < 1e-14);
assert(norm(se-sqrt(diag(J*V*J'))) < 1e-14 && p(3) > .999);
[~,clippedSE] = hybrid_cutoffs([-2;20],V(1:2,1:2),10);
assert(clippedSE(2) == sqrt(V(1,1)));
for fullCov = 0:2
    e = struct('NVarA',3,'NLatent',2,'NVarM',2,'FullCov',fullCov);
    expected = [12 15 22];
    assert(isequal(hybrid_lv_indices('HMXL',e),expected(fullCov+1)+(1:6)));
    assert(isequal(hybrid_lv_indices('HMNL',e),9+(1:6)));
end
fprintf('HYBRID_OPTIMIZATION_UNIT_PASS\n');
    function [f,g,H] = quadratic(b)
        v = b-target; f = .5*v'*Q*v; g = Q*v; H = Q;
    end
    function stop = checkCallback(b,values,~)
        assert(numel(b) == 3 && b(3) == 7);
        if isfield(values,'gradient') && ~isempty(values.gradient), assert(numel(values.gradient) == 3); end
        seen = true; stop = false;
    end
end
