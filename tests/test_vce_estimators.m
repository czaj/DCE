function test_vce_estimators(outDir)
if ~isfolder(outDir), mkdir(outDir); end
settings = parallel.Settings; previous = settings.Pool.AutoCreate;
restore = onCleanup(@() setAuto(settings,previous)); %#ok<NASGU>
settings.Pool.AutoCreate = false; assert(isempty(gcp('nocreate')));
for model = {'MDCEV','MMDCEV'}
    name = model{1}; c = vce_fixture(name,3,1,1,1,8,16,3,3); e = c.opt;
    I = c.data; I.Y = I.Y(:); I.priceMat = I.priceMat(:);
    if strcmp(name,'MMDCEV')
        I.Xa = reshape(permute(I.Xa,[1 3 2]),e.NAlt*e.NCT*e.NP,e.NVarA);
        I.Xm = repelem(I.Xm',e.NAlt*e.NCT,1); I = rmfield(I,'err');
    else
        I.Xm = repelem(I.Xm,e.NAlt,1);
    end
    I.Xu = repelem(I.Xu(:,2:end),e.NAlt,1);
    e.Display = 0; e.ConstVarActive = 0; e.ApproxHess = 1; e.NSdSim = 100;
    e.Seed1 = 123; e.Seed2 = 456; e.Draws = 1; e.NObs = e.NCT*e.NP;
    e.ProjectName = ['test_',name]; e.SaveTxtOutput = 0; e.OutputDir = outDir;
    e.BActive(end) = 0; old.(name) = struct('b0',c.b);
    oo = optimset('Algorithm','quasi-newton','GradObj','on','Hessian','off',...
        'MaxIter',1,'Display','off','FinDiffType','central');
    for fix = [1 3 5]
        e.HessEstFix = fix; r = feval(name,I,old,e,oo);
        assert(r.EstimOpt.NumGrad == 0 && all(isfinite(r.hess),'all'));
        input = r.INPUT;
        if strcmp(name,'MMDCEV'), input.err = generateRandomDraws(r.EstimOpt)'; end
        fun = @(b) feval(['probs_',lower(name)],input,r.EstimOpt,b);
        W = input.W; if strcmp(name,'MMDCEV'), W = W(1:e.NCT:end); end
        assert(isequal(input.W,I.W),'Weights changed.');
        active = r.EstimOpt.BActive ~= 0;
        if fix == 1
            [~,j] = fun(r.bhat); j = j(:,active).*W; expected = j'*j;
        else
            expected = hessianFromGradient(@(b) gradient(fun,W,b),r.bhat,r.EstimOpt.BActive);
            expected = expected(active,active);
        end
        error = max(abs(r.hess-expected),[],'all')/max(1,max(abs(expected),[],'all'));
        assert(error < 5e-4,'Public Hessian mismatch %s fix %d: %.3g',name,fix,error);
        save(fullfile(outDir,sprintf('%s_hessian%d.mat',name,fix)),'r','expected','error');
        fprintf('VCE_PUBLIC_PASSED %s HessEstFix=%d error=%.3g\n',name,fix,error);
    end
end
end
function g = gradient(fun,W,b)
[~,j] = fun(b); g = sum(j.*W,1)';
end
function setAuto(settings,value)
settings.Pool.AutoCreate = value;
end
