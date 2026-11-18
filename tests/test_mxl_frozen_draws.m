function comparison = test_mxl_frozen_draws(outDir)
% Public MXL -> frozen refit, compared with a separately coded binary panel integral.
% Request an output to also time the correlated WTP example at two draw budgets.
if nargin<1, outDir = tempname; end
if ~isfolder(outDir), mkdir(outDir); end
for wtp = 0:1
    for fullCov = 0:1
        rng(918367);
        N = 180; T = 12; d = 2;
        X = randn(2,T,N,d);
        X(:,:,:,2) = -2*rand(2,T,N);
        rho = .65*fullCov;
        logSD = .65+1.75*wtp;
        factor = [1.1 0;logSD*rho logSD*sqrt(1-rho^2)];
        latent = [0.8;-.3]+factor*randn(d,N);
        money = exp(latent(2,:));
        taste = latent(1,:).*money.^wtp;
        V = X(:,:,:,1).*reshape(taste,1,1,N)+X(:,:,:,2).*reshape(money,1,1,N);
        [~,choice] = max(V-log(-log(rand(2,T,N))),[],1);
        Y = cat(1,choice==1,choice==2);
        I = struct('Y',double(Y(:)),'Xa',reshape(X,[],d));
        observed = true(T,N);
        if ~wtp && ~fullCov
            observed(1,1:3:end) = false;
            I.Y(repelem(~observed(:),2)) = NaN;
        end
        E = struct('NP',N,'NCT',T,'NAlt',2,'NamesA',{{'attribute';'negative cost'}}, ...
            'Display',0,'OutputDir',outDir,'SaveTxtOutput',0,'ProjectName','frozen');
        [I,old,E,O] = DataCleanDCE(I,E);
        E.WTP_space = wtp; E.FullCov = fullCov; E.Dist = [0 1];
        E.RobustStd = 0; E.HessEstFix = 5; E.EffectiveMoments = 0; E.Drawskeep = 1;
        E.Scores = 0; E.GPU = 'cpu'; E.Seed1 = 17431; E.NRep = 512;
        E.ImportanceSampling = 'fixed'; E.Display = wtp && fullCov;
        E.OutputFilePrefix = 'frozen_full';
        O.Algorithm = 'quasi-newton'; O.Hessian = 'off'; O.GradObj = 'on'; O.Display = 'off';
        start = [.4;0;-.9;.8]; field = 'MXL_d';
        if fullCov, start = [.4;0;-.9;-.1;.8]; field = 'MXL'; end
        old.(field).b0 = start;
        if ~wtp, I.W = 0.8 + (1:N)'/(2*N); end
        result = MXL(I,old,E,O);
        fit = result.IS;
        assert(fit.success && strcmp(fit.Method,'fixed') && numel(fit.Rounds)==1 && result.EstimOpt.Draws==6);
        assert(abs(fit.StartLL-result.LL)<1e-7 && fit.LL>=fit.StartLL-1e-7);
        assert(norm(sum(fit.jacobian,1)'-fit.g)<1e-8 && fit.diagnostics.MaxNewtonStepSE<.01);
        assert(all(fit.ESS>=1-1e-8 & fit.ESS<=fit.NRep+1e-8) && all(isfinite(fit.RefreshStepSE)));
        R = fit.NRep;
        mask = logical(eye(d)); if fullCov, mask = tril(true(d)); end
        pilot = result.bhat;
        L = zeros(d); L(mask) = pilot(d+1:end);
        Z = reshape(pilot(1:d)+L*result.err,d,R,N);
        % Direct independent calculation with determinants and weighted per-person log integrals.
        s = reshape(2*double(choice==1)-1,T,N);
        dx = reshape(X(1,:,:,:)-X(2,:,:,:),T,N,d).*s;
        logP = zeros(R,N);
        for n = 1:N
            beta = Z(:,:,n);
            beta(2,:) = exp(beta(2,:)); beta(1,:) = beta(1,:).*beta(2,:).^wtp;
            v = -reshape(dx(:,n,:),T,d)*beta;
            v = v(observed(:,n),:);
            logP(:,n) = -sum(max(v,0)+log1p(exp(-abs(v))),1)';
        end
        base = logP-reshape(logNormal(reshape(Z,d,[]),pilot(1:d),L),R,N);
        criterion = @(b) integralAt(b,Z,base,mask,I.W);
        [value,panels] = criterion(fit.bhat);
        assert(abs(value+fit.LL)<1e-7 && norm(panels-fit.LLdetailed,Inf)<1e-8);
        k = numel(fit.bhat); step = 8e-5*(1+abs(fit.bhat)); A = diag(step);
        H = zeros(k);
        for i = 1:k
            up = fit.bhat+A(:,i); down = fit.bhat-A(:,i);
            [~,plus] = criterion(up); [~,minus] = criterion(down);
            assert(norm((plus-minus)/(2*step(i))-fit.jacobian(:,i),Inf)<3e-5);
            for j = 1:i
                h = (criterion(up+A(:,j))-criterion(up-A(:,j))-criterion(down+A(:,j))+criterion(down-A(:,j)))/(4*step(i)*step(j));
                H(i,j) = h; H(j,i) = h;
            end
        end
        se = sqrt(diag(H\eye(k)));
        assert(max(abs(se./fit.std-1))<1e-3);
        if fullCov
            L(mask) = fit.bhat(d+1:end);
            assert(norm(L*L'-fit.DetailsVcov,'fro')<1e-10 && abs(fit.DetailsVcor(1,2)-L(2,1)/norm(L(2,:)))<1e-10);
            rho = @(b) b(4)/hypot(b(4),b(5));
            delta = zeros(k,1);
            for j = 1:k, delta(j) = (rho(fit.bhat+A(:,j))-rho(fit.bhat-A(:,j)))/(2*step(j)); end
            assert(abs(fit.LatentCorrelations(1,4)/sqrt(delta'*(H\delta))-1)<1e-3);
            assert(abs(fit.LatentCorrelations(1,3)-rho(fit.bhat))<1e-12 && fit.LatentCorrelations(1,5)<=1);
            assert(norm(fit.DetailsV(:,1)-sqrt(diag(L*L')))<1e-10 && all(isnan(fit.DetailsV(:,4))));
        end
        if E.Display
            sheet = readcell(fit.OutputFile,'Sheet','Results');
            assert(abs(sheet{4,5}-fit.std(1))<1e-9 && any(strcmp(sheet(:,1),'Fresh IS draws')));
            row = find(strcmp(sheet(:,1),'negative cost / attribute'),1);
            assert(~isempty(row) && abs(sheet{row,3}-fit.LatentCorrelations(1,3))<1e-9);
            assert(abs(sheet{row,5}-fit.LatentCorrelations(1,4))<1e-9);
        end
        fprintf('FROZEN_PASS WTP=%d full=%d LL=%.8f SE discrepancy=%.4g fresh=%.3f SE\n',wtp,fullCov,fit.LL,max(abs(se./fit.std-1)),max(fit.RefreshStepSE));
    end
end
if nargout
    comparison = compareCorrelated(I,old,E,O);
end
% New budget and fresh scramble are reproducible independently of the pilot grid.
E.Display = 0; E.ISNRep = 1024; E.ISSeed = 8765;
b = MXL(I,old,E,O);
assert(b.IS.success && b.IS.NRep==1024 && b.IS.Seed==8765 && abs(b.IS.StartLL-b.LL)>1e-7);
assert(abs(b.LL-result.LL)<1e-8 && max(abs(b.IS.std-fit.std))>1e-6);
% A coefficient with no information remains free; it cannot yield valid inference on a reduced model.
E.Display=0; E.WTP_space=0; E.FullCov=0; E.ImportanceSampling=0;
I.Xa(:,1)=0; old.MXL_d.b0=[.4;0;.9;.8];
flat = MXL(I,old,E,O);
assert(~flat.HessDiagnostics.ValidInference && all(isfinite(flat.bhat)) && all(isnan(flat.std)));
assert(all(flat.EstimOpt.BActive) && flat.diagnostics.NNonFiniteSE==4 && ~isfield(flat,'IS'));
fprintf('ALL_MXL_FROZEN_TESTS_PASS\n');
end

function comparison = compareCorrelated(I,old,E,O)
% Time complete public calls after warm-up; the fixed strategies include the 512-draw pilot.
E.Display = 0; E.Drawskeep = 0;
strategy = ["SML512","fixed512","SML4096","fixed4096"];
records = cell(24,1);
row = 0;
for seed = [17431 71629]
    E.Seed1 = seed;
    for rep = 1:3
        for index = circshift(1:4,rep-1)
            fixed = startsWith(strategy(index),"fixed");
            budget = 512; if index>2, budget=4096; end
            E.NRep = budget; E.ImportanceSampling = 0;
            if fixed
                E.NRep = 512; E.ISNRep = budget; E.ImportanceSampling = 'fixed';
            end
            clock = tic;
            pilot = MXL(I,old,E,O);
            elapsed = toc(clock);
            fit = pilot; extra = 0;
            if fixed
                assert(pilot.IS.success && strcmp(pilot.IS.Method,'fixed'));
                fit = pilot.IS; extra = fit.EstimationSeconds;
            end
            assert(fit.exitf>0 && fit.diagnostics.MaxNewtonStepSE<.01 && all(isfinite(fit.std)));
            pair = fit.LatentCorrelations(1,:);
            row = row+1;
            records{row} = struct('Seed',seed,'Repeat',rep,'Method',strategy(index), ...
                'Seconds',elapsed,'RefitSeconds',extra,'LL',fit.LL, ...
                'MeanWTP',fit.bhat(1),'SEMeanWTP',fit.std(1),'SDWTP',fit.DetailsV(1,1), ...
                'SESDWTP',fit.DetailsV(1,3),'Rho',pair(3),'SERho',pair(4), ...
                'MeanLogScale',fit.bhat(2),'SEMeanLogScale',fit.std(2), ...
                'SDLogScale',fit.DetailsV(2,1),'SESDLogScale',fit.DetailsV(2,3));
        end
    end
end
comparison = struct2table(vertcat(records{:}));
% Independent 224-by-224 Gauss-Hermite refit on this dataset (see DRAWS.md).
reference = [0.17269118,0.14053584,0.15751111,0.24329305,0.35183601];
fixed = comparison.Method == "fixed4096";
se = comparison{fixed,{'SEMeanWTP','SESDWTP','SERho','SEMeanLogScale','SESDLogScale'}};
assert(max(abs(se./reference-1),[],'all')<0.10,'High-budget fixed refit departed from quadrature.');
end

function v = logNormal(points,mu,L)
cov = L*L';
delta = points-mu;
v = -(sum(delta.*(cov\delta),1)+log(det(cov))+numel(mu)*log(2*pi))/2;
end

function [value,panels] = integralAt(b,Z,base,mask,W)
[d,R,N] = size(Z);
L = zeros(d); L(mask) = b(d+1:end);
logW = base+reshape(logNormal(reshape(Z,d,[]),b(1:d),L),R,N);
m = max(logW,[],1);
panels = -W(:).*(m+log(mean(exp(logW-m),1)))';
value = sum(panels);
end
