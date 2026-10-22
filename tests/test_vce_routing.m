function test_vce_routing
for model = {'MDCEV','MMDCEV'}
    c = vce_fixture(model{1},3,1,1,1);
    for mode = {'cpu','auto','gpu'}
        c.opt.GPU = mode{1};
        for change = 0:2
            if change == 1, c.b(1) = c.b(1)+.1; end
            if change == 2, c.data.Xa(1) = c.data.Xa(1)+.2; end
            d = vce_prepare(c.model,c.data,c.opt);
            [fc,gc] = vce_evaluate(d,c.b,2,true,0,[]);
            for repeat = 1:2
                [f,g,backend] = vce_likelihood(c.model,c.data,c.opt,c.b,true);
                assert(max(abs(f-fc)) < 1e-8 && max(abs(g-gc),[],'all') < 1e-6);
            end
            fprintf('VCE routing %s %s -> %s change=%d\n',c.model,mode{1},backend,change);
        end
    end
    c.opt.GPU = 'cpu'; c.opt.BActive(end) = 0;
    oo = struct('GradObj','on','Hessian','user-supplied','FinDiffType','central');
    wrapper = ['LL_',lower(c.model),'_MATlike'];
    c.opt.NumGrad = 0; [fa,ga,ha] = feval(wrapper,c.data,c.opt,oo,c.b);
    c.opt.NumGrad = 1; [fn,gn,hn] = feval(wrapper,c.data,c.opt,oo,c.b);
    assert(abs(fa-fn) < 1e-10 && max(abs(ga-gn)) < 1e-6 && max(abs(ha-hn),[],'all') < 1e-5);
    oo.GradObj = 'off'; fv = feval(wrapper,c.data,c.opt,oo,c.b); assert(abs(fv-fa) < 1e-10);
    if strcmp(c.model,'MMDCEV')
        c.data.W(2) = c.data.W(1)+1; failed = false;
        try, feval(wrapper,c.data,c.opt,oo,c.b); catch ex, failed = strcmp(ex.identifier,'MMDCEV:Weights'); end
        assert(failed,'Inconsistent panel weights were accepted.');
    end
end
for bad = {'invalid',struct('unexpected',1)}
    c.opt.GPU = bad{1}; failed = false;
    try
        vce_likelihood(c.model,c.data,c.opt,c.b,true);
    catch ex
        failed = strcmp(ex.identifier,'DCE:VCE:GPUOption');
    end
    assert(failed,'Invalid GPU option was accepted.');
end
c.opt.GPU = 'cpu'; c.opt.VCEBlockSize = 0; failed = false;
try
    vce_likelihood(c.model,c.data,c.opt,c.b,true);
catch ex
    failed = strcmp(ex.identifier,'DCE:VCE:BlockSize');
end
assert(failed,'Invalid block size was accepted.');
fprintf('VCE_ROUTING_PASSED\n');
end
