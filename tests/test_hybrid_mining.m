function report = test_hybrid_mining(fitFile,outDir)
if ~isfolder(outDir), mkdir(outDir); end
s = load(fitFile,'fit'); fit = s.fit; I = fit.INPUT; o = fit.EstimOpt;
o.NRep = 64; o.GPU = 'cpu';
rng(o.Seed1); h = sobolset(o.NVarA+o.NLatent,'Skip',o.HaltonSkip,'Leap',o.HaltonLeap);
h = scramble(h,'MatousekAffineOwen'); err = norminv(net(h,o.NP*o.NRep))';
args = {I.YY,I.XXa,I.Xm,I.Xs,I.Xstr,I.Xmea,I.Xmea_exp,err,o,fit.bhat};
settings = parallel.Settings; previousAuto = settings.Pool.AutoCreate;
restoreAuto = onCleanup(@() setAutoCreate(settings,previousAuto)); %#ok<NASGU>
settings.Pool.AutoCreate = false;
assert(isempty(gcp('nocreate')));
[f,g] = LL_hmxl(args{:});
old = args; old{9}.HybridBlockSize = 0; [fo,go] = LL_hmxl(old{:});
fprintf('Mining baseline errors: value %.12g, gradient %.12g\n',relative(f,fo),relative(g,go));
save(fullfile(outDir,'comparison.mat'),'f','g','fo','go','args');
unclamped = f < -log(o.ProbMin)-1e-8;
assert(relative(f,fo) < 1e-8 && relative(g(unclamped,:),go(unclamped,:)) < 1e-6,'Mining baseline mismatch.');
assert(all(g(~unclamped,:) == 0,'all'),'A floored objective has zero derivative.');
J = zeros(size(g));
for j = 1:numel(fit.bhat)
    step = 1e-5*max(1,abs(fit.bhat(j))); ap = args; am = args;
    ap{10}(j) = ap{10}(j)+step; am{10}(j) = am{10}(j)-step;
    J(:,j) = (LL_hmxl(ap{:})-LL_hmxl(am{:}))/(2*step);
end
assert(relative(g,J) < 5e-6,'Mining finite difference mismatch.');
active = o.BActive; active(end) = 0;
H = hessianFromGradient(@gradient,fit.bhat,active);
assert(all(isfinite(H),'all') && max(abs(H-H'),[],'all') < 1e-10);
assert(all(H(~active,:) == 0,'all') && all(H(:,~active) == 0,'all'));
% The optimizer wrapper must keep weights and fixed-parameter scores intact.
oo = struct('GradObj','on','Hessian','user-supplied','FinDiffType','central');
o.BActive = active;
[objective,score,bhhh] = LL_hmxl_MATlike(args{1:8},I.W,o,oo,fit.bhat);
jw = g.*I.W; jw(:,~active) = 0;
assert(abs(objective-sum(f.*I.W)) < 1e-8 && relative(score,sum(jw,1)') < 1e-8);
assert(relative(bhhh,jw'*jw) < 1e-8);
report = struct('NP',o.NP,'NRep',o.NRep,'Parameters',numel(fit.bhat),...
    'ValueError',relative(f,fo),'GradientError',relative(g(unclamped,:),go(unclamped,:)),...
    'FlooredRespondents',sum(~unclamped),'LegacyFloorGradient',max(abs(go(~unclamped,:)),[],'all'),...
    'FiniteDifferenceError',relative(g,J));
save(fullfile(outDir,'mining.mat'),'report','f','g','J','H');
fid = fopen(fullfile(outDir,'mining.json'),'w'); fprintf(fid,'%s',jsonencode(report)); fclose(fid);
disp(report); fprintf('MINING_TEST_PASSED\n');
    function value = gradient(b)
        a = args; a{end} = b; [~,scores] = LL_hmxl(a{:}); value = sum(scores.*I.W,1)';
    end
end

function r = relative(a,b)
r = max(abs(a-b),[],'all')/max(1,max(abs(b),[],'all'));
end

function setAutoCreate(settings,value)
settings.Pool.AutoCreate = value;
end
