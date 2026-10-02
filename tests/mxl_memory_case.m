function C = mxl_memory_case(name,replicationRoot)
% Shared input for correctness checks and allocation benchmarks.
% Names: demo_pref_diag, demo_pref_full, demo_wtp_diag, demo_wtp_full, CH, pooled.
repo = fileparts(fileparts(mfilename('fullpath')));
if nargin < 2 || isempty(replicationRoot)
    replicationRoot = fullfile(fileparts(fileparts(fileparts(repo))),...
        'lasy','replication_package');
end
addpath(fullfile(repo,'MXL'),fullfile(repo,'MNL'),...
    genpath(fullfile(fileparts(repo),'tools')));
name = char(name);
C.Name = name;
C.Published = [];
E = struct('SaveTxtOutput',0,'OutputDir',tempdir,'CheckSeparation',0);

if startsWith(name,'demo_')
    assert(ismember(name,{'demo_pref_diag','demo_pref_full',...
        'demo_wtp_diag','demo_wtp_full'}),'Unknown demo case: %s',name);
    dataFile = fullfile(repo,'DCE_demo','NEWFOREX_DCE_demo.mat');
    D = load(dataFile);
    if isfield(D,'Choice')
        INPUT.Y = D.Choice(:);
    else
        INPUT.Y = D.Y(:);
    end
    assert(numel(INPUT.Y) == numel(D.SQ),'Demo Y must be a long-format indicator.');
    if isfield(D,'SKIP'), INPUT.MissingInd = D.SKIP(:); end
    INPUT.Xa = [D.SQ,D.GOS,D.CEN,D.VIS == 2,D.VIS == 1,-D.FEE/4/10];
    E.NamesA = {'Status quo';'GOS';'CEN';'VIS=2';'VIS=1';'-Cost (10 EUR)'};
    E.NCT = 12;
    E.NAlt = 3;
    E.NP = numel(INPUT.Y)/(E.NCT*E.NAlt);
    E.Dist = [0 0 0 0 0 1];
    E.WTP_space = double(contains(name,'_wtp_'));
    E.FullCov = double(endsWith(name,'_full'));
    E.Display = 0;
    [INPUT,C.Results,E,C.OptimOpt] = DataCleanDCE(INPUT,E);
    if E.FullCov == 0
        C.b = [-1.0862014791638248;.6423214063535723;.8650703151893991;...
            .2607394211805801;.1073491359204689;-.0992110619362219;...
            3.0447049702322388;.8722468228131102;.8938160900831796;...
            .6118102989399823;-.2365315344749905;1.2642393272248089];
    else
        C.b = [-1.1972997845872539;.5646621201319871;.8406514728208795;...
            .0564480918844042;-.0420225209056210;-.2495359314908631;...
            3.1494409625688204;-.5542040007017031;-.4605229352346021;...
            -.5237055639407464;-.3886718922677192;-.4440679268877266;...
            .8991462416618976;.4315092142776225;-.0979503116131898;...
            -.2494899587824761;.1923782446226067;.9082669077418934;...
            .2336754749973365;.0949517602128645;.0474132249817506;...
            .7253124386918487;.3323448145561994;-.3466234016364159;...
            -.2353464620237843;.3358189295279997;.8981782979719334];
    end
else
    assert(ismember(name,{'CH','pooled'}),'Unknown forest case: %s',name);
    addpath(fullfile(replicationRoot,'code'));
    DATA = load_forest_data(fullfile(replicationRoot,'data','forest_dce_data.mat'));
    S = load(fullfile(replicationRoot,'results','published','estimates_published.mat'),'PUB');
    filter = DATA.estim_sample;
    if strcmp(name,'CH')
        filter = filter & DATA.country_model == "CH";
        C.Published = S.PUB.CH.MXL;
        expectedNP = 644;
    else
        C.Published = S.PUB.ALL_paper.MXL;
        expectedNP = 8940;
    end
    INPUT.Y = DATA.Y(filter);
    INPUT.Xa = DATA.Xa(filter,:);
    INPUT.W = DATA.Weq(filter,:);
    E.NamesA = DATA.NamesA;
    E.NCT = 12;
    E.NAlt = 3;
    E.NP = numel(INPUT.Y)/(E.NCT*E.NAlt);
    E.Display = 0;
    [INPUT,C.Results,E,C.OptimOpt] = DataCleanDCE(INPUT,E);
    assert(E.NP == expectedNP,'Unexpected estimation sample: %s',name);
    E.FullCov = 1;
    E.WTP_space = 1;
    E.NRep = 1000;
    E.Dist = [zeros(1,16),1];
    C.b = C.Published.bhat(:);
end

% These are the same indices and shapes prepared by MXL before LL_mxl.
E.NVarA = size(INPUT.Xa,2);
if ~isfield(INPUT,'Xm'), INPUT.Xm = zeros(numel(INPUT.Y),0); end
if ~isfield(INPUT,'Xs'), INPUT.Xs = zeros(numel(INPUT.Y),0); end
E.NVarM = size(INPUT.Xm,2);
E.NVarS = size(INPUT.Xs,2);
E.NVarNLT = 0;
E.NLTVariables = [];
E.NLTType = [];
E.Johnson = 0;
E.Triang = [];
E.ExpB = [];
E.NumGrad = 0;
E.Dist = E.Dist(:)';
E.WTP_matrix = [];
if E.WTP_space > 0, E.WTP_matrix = E.NVarA*ones(1,E.NVarA-1); end
C.YY = reshape(INPUT.Y,E.NAlt*E.NCT,E.NP);
C.XXa = permute(reshape(INPUT.Xa,E.NAlt*E.NCT,E.NP,E.NVarA),[1 3 2]);
C.XXm = reshape(INPUT.Xm',[E.NVarM,E.NAlt*E.NCT,E.NP]);
E.mCT = false;
for n = 1:E.NP
    for k = 1:E.NVarM
        x = C.XXm(k,:,n);
        x = x(isfinite(x));
        E.mCT = E.mCT || (~isempty(x) && any(x ~= x(1)));
    end
end
if E.mCT
    C.XXm = INPUT.Xm';
else
    C.XXm = reshape(C.XXm(:,1,:),E.NVarM,E.NP);
end
C.Xs = INPUT.Xs;
C.W = INPUT.W(:);
index = tril(ones(E.NVarA));
index(index == 1) = 1:sum(1:E.NVarA);
E.DiagIndex = diag(index);
E.indx1 = [];
E.indx2 = [];
for k = 1:E.NVarA
    E.indx1 = [E.indx1,k:E.NVarA];
    E.indx2 = [E.indx2,k*ones(1,E.NVarA+1-k)];
end
C.EstimOpt = E;
C.INPUT = INPUT;

oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng));
rng(E.Seed1);
assert(E.Draws == 6,'These fixtures use the scrambled Sobol specification.');
sequence = sobolset(E.NVarA,'Skip',E.HaltonSkip,'Leap',E.HaltonLeap);
sequence = scramble(sequence,'MatousekAffineOwen');
draws = icdf('Normal',net(sequence,E.NP*E.NRep),0,1);
draws(:,E.Dist == -1) = 0;
C.err = draws';
end
