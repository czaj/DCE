
% DCE demo script - 2017-03-06
% (CC BY 4.0) czaj.org

% The data used here was used for, e.g.:
% Czajkowski, M., Barczak, A., Budzinski, W., Giergiczny, M., and Hanley, N., 2016. Preference and WTP stability for public forest management. Forest Policy and Economics, 71:11-22.

% The script requires:
% https://github.com/czaj/dce
% https://github.com/czaj/tools
% (they need to be downloaded and included in Matlab paths (with subfolders))

clear all
clc

global B_backup; % this is for storing B (parameters) vector in case iterations are interrupted with ctrl-c


%% ****************************  loading data  ****************************


EstimOpt.DataFile = ('NEWFOREX_DCE_demo.mat'); % specify the name of the data file
DATA = load(EstimOpt.DataFile); % load data, store it in a structure DATA


%% ****************************  data transformations  ****************************



%% ****************************  model specification  ****************************


DATA.Y = DATA.Y; % exlpanained variable (choice indicator)

DATA.Xa = [DATA.SQ, DATA.GOS, DATA.CEN, DATA.VIS==2, DATA.VIS==1, -DATA.FEE/4/10]; % attributes
EstimOpt.NamesA = {'Status quo';'Passive protection of all commertial forests';'Passive protection of all natural regeneration forests';'Reducing the number of visitors to 7500 per day';'Reducing the number of visitors to 5000 per day';'-Cost (10 EUR)'}; % Names of the attributes

% DATA.Xm = []; 
% EstimOpt.NamesM = {''}; 

% DATA.Xs = [];  
% EstimOpt.NamesS = {''};

% INPUT.Xt = []; 
% EstimOpt.NamesT = {''}; 

% DATA.Xc = []; 
% EstimOpt.NamesC = {''}; 

% DATA.Xmea = []; % 
% EstimOpt.NamesMea = {''};

% DATA.Xmea_exp = [];
% EstimOpt.NamesMea_exp = {''};

% DATA.Xstr = [];
% EstimOpt.NamesStr = {''};

EstimOpt.ProjectName = 'demo'; % name of the project used for naming xls results files


%% ****************************  specifying input ****************************


DATA.filter = ones(size(DATA.Y)) == 1; % modify to include only some observations, e.g. model for male respondents only

INPUT.Y = DATA.Y(DATA.filter); 
INPUT.Xa = DATA.Xa(DATA.filter,:);
% INPUT.Xm = DATA.Xm(DATA.filter,:); % Explanatory variables of random parameters means / inreactions
% INPUT.Xs = DATA.Xs(DATA.filter,:); % Explanatory variables of scale
% INPUT.Xt = DATA.Xt(DATA.filter,:); % Explanatory variables of scale variance (GMXL model)
% INPUT.Xc = DATA.Xc(DATA.filter,:); % Explanatory variables of class membership (Latent class models)
% INPUT.Xmea = DATA.Xmea(DATA.filter,:); % Measurment equations variables (Hybrid models)
% INPUT.Xmea_exp = DATA.Xmea_exp(DATA.filter,:); % Additional covariates explaining mesurment equations (Hybrid models)
% INPUT.Xstr = DATA.Xstr(DATA.filter,:); % Structural equations variables (Hybrid models)

% INPUT.MissingInd = DATA.SKIP(DATA.filter,:); % use to indicate missing observations, e.g. choice tasks with no answer


%% ****************************  sample characteristics ****************************


EstimOpt.NCT = 12; % Number of choice tasks per person 
EstimOpt.NAlt = 3; % Number of alternatives
EstimOpt.NP = length(INPUT.Y)/EstimOpt.NCT/EstimOpt.NAlt; % 789; % Number of respondents


%% **************************** estimation and optimization options ****************************


% EstimOpt.eps = 1.e-9; % overall precision level

[INPUT, Results, EstimOpt, OptimOpt] = DataCleanDCE(INPUT,EstimOpt);

% EstimOpt.NRep = 1e4; % number of draws for numerical simulation
% OptimOpt.MaxIter = 1e3; % maximum number of iterations
% OptimOpt.Algorithm = 'trust-region'; % 'quasi-newton'
% EstimOpt.NumGrad = 1; % 0
% OptimOpt.GradObj = 'off'; % 'off'
% OptimOpt.FinDiffType = 'central'; % 'forward'
% OptimOpt.Hessian = 'user-supplied'; % 'off'
% EstimOpt.HessEstFix = 1; % 0 = use optimization Hessian, 1 = use jacobian-based (BHHH) Hessian, 2 - use high-precision jacobian-based (BHHH) Hessian 3 - use numerical Hessian
% EstimOpt.ApproxHess = 0;

% Estimopt.RealMin = 1; % use in the case of numerical errors (possibly for a few iterations only)


%% ****************************     MNL     ****************************


% EstimOpt.WTP_space = 1; % number of monetary attributes for WTP-space estimation (need to come last in Xa)
% EstimOpt.WTP_matrix = []; % specify which monetary parameter is used for which non-monetary attribute for WTP-space models
% EstimOpt.NLTType = 1 % 1 = Box-Cox transformations; 2 = Yeo-Johnson transofmarions (MNL and MXL only)
% EstimOpt.NLTVariables = []; % choose Xa to be non-linearly transformed (MNL and MXL only)

Results.MNL.b0 = [0.159954254207451;0.486588294026178;0.609085577346533;0.188275582287589;0.0682110004987603;0.569893495060325]; % provide starting values for the MNL model

Results.MNL = MNL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     MXL_d     ****************************


% EstimOpt.FullCov = 0; % Specify if random parametrs are correlated or not (0 = default)
% EstimOpt.Dist = zeros(size(INPUT.Xa,2),1); % choose Xa parameters' distributions (0 - normal, -1 - constant, 1 - lognormal, 2 - Spike, 3 - triangular, 4 - Weibull, 5 - sinh-arcsinh, 6 - Johnson Sb, 7 - Johnson Su); 
% EstimOpt.Scores = 1; % Predict individual-specific paramter / WTP scores using bayes rule
 
Results.MXL_d.b0 = [-0.740090971865778;0.637808759156240;0.862255786852351;0.259514488196178;0.119545771432164;1.20126366401120;3.26204813109697;0.889290371447226;0.921208308805706;0.585533284920086;0.250373926072790;1.03957224553947]; % provide starting values for the MXL_d model

Results.MXL_d = MXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     MXL     ****************************


EstimOpt.FullCov = 1;

Results.MXL.b0 = [-0.990816054697052;0.510612251809615;0.790917452248288;-0.0263706789518912;-0.0521799755969684;1.22507452590470;2.90674179266798;-0.487180784239663;-0.414669765632886;-0.558920458472267;-0.383750973613273;-0.253255599241274;0.879987387494616;0.397422267248830;-0.146603943803091;-0.266287907573810;0.0817081768967246;0.921880514420358;0.226640541893134;0.123818012597966;0.0121953367738178;0.615768549750211;0.347581885106293;-0.682792855103666;0.150392172041538;-0.0995369723675956;0.785555928723963]; % provide starting values for the MXL model

Results.MXL = MXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     GMXL_d     ****************************


EstimOpt.FullCov = 0;
% EstimOpt.Gamma0 = 0.5; % 0 and 1 impose GMXL type II and I, respectively. Value in between allows for gamma to be freely estimated.

Results.GMXL_d.b0 = [-1.06455191723154;0.915763550037752;1.25973862236508;0.421581388848477;0.238181305805212;1.77569286454493;4.44099589447681;1.05012841199027;0.975206824398307;0.723397561577961;-0.145911717636068;1.48303763037042;0.923169871577123;0]; % provide starting values for the GMXL_d model

Results.GMXL_d = GMXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     GMXL     ****************************


EstimOpt.FullCov = 1;

Results.GMXL.b0 = [-1.32448971675538;0.650407429819206;1.15273787184460;-0.0739900072709748;-0.0853291917786047;1.75908163393569;3.64076267564284;-0.620949814952573;-0.367734713782486;-0.658710526284171;-0.506934789113255;-0.123254247784516;1.03754456555579;0.498440815134805;-0.222482192403819;-0.330764378530705;0.0177926073795942;1.07316872059076;0.219976825369281;0.117670320704554;0.227275431831907;0.758181720965921;0.392905084260305;-0.997554368267189;0.0724136026286958;-0.577421700883042;0.771441424386692;0.802607064191922;0]; % provide starting values for the GMXL model

Results.GMXL = GMXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     LC     ****************************


% EstimOpt.NClass = 2; % number of latent classes
% EstimOpt.ClassScores = 1; % Predict individual-specific class membership probabilities using Bayes rule
% EstimOpt.BActiveClass = ones(EstimOpt.NVarA,1); % use 0 to impose coefficient equality for all latent classes

Results.LC.b0 = [0.910707942150065;0.304347621653407;0.0275411336605089;-0.463670945583358;-0.173647468159989;3.35386438616489;-1.43265608444431;0.385891077706399;0.634607481281206;0.400736597552428;0.326055705149880;0.438582404864380;-0.425398435055134];

Results.LC = LC(INPUT,Results,EstimOpt,OptimOpt);

xlswrite('LC_results_demo.xls', Results.LC.R_out)


%% ****************************     LCMXL     ****************************


EstimOpt.FullCov = 0;
% EstimOpt.Scores = 1; % Predict individual-specific paramters / class membership probabilities using Bayes rule

Results.LCMXL_d.b0 = [0.355247459379893;0.284008208543304;0.363107885792539;-0.100846141867809;-0.152063399154145;1.02066359651063;-3.68013741701911;1.41381415235096;2.05731373376956;1.11161663077177;0.705747412459405;1.80969785093245;2.44471544846888;0.453259807660528;-0.340508482472775;0.336650542081059;0.322439038465852;0.858761681902807;1.73271002124883;1.31859132566682;1.38848708832864;-0.875152589043110;-0.0600842052577111;1.39048282621374;0.625923593008172];

Results.LCMXL_d = LCMXL(INPUT,Results,EstimOpt,OptimOpt);

xlswrite('LCMXL_d_results_demo.xls', Results.LCMXL_d.R_out)


%% ****************************     LCMXL     ****************************


EstimOpt.FullCov = 1;

Results.LCMXL.b0 = [0.278937432867506;0.291685483501190;0.433222289479560;-0.188629458329337;-0.211552385198981;0.961867184080436;-3.91199316818707;1.54484147484562;2.51646857877832;0.661558953759638;0.530048394522240;2.09872944791070;2.38295904299256;-0.189652838261572;0.0893664139057544;-0.242452835824827;-0.0988226206906153;0.0787048588021720;0.510419009609970;0.429373135410626;0.0984778037994452;-0.207044139857037;0.267041687173508;0.0528162367313782;-0.661637708068180;-0.440408457308465;0.486718693283753;0.0341381103896203;-0.0480276684030333;-0.279327145876270;-0.299693403252834;0.201197015075489;0.487351572930574;2.15831632772265;-0.171233942614560;0.720082525467859;-0.588052507283012;-0.356246077020086;-0.479011542123428;1.74598170056253;0.417345679713036;-0.267628531774601;-0.366892405168353;-0.430206693338263;1.08724269832352;0.824320380698708;0.459788099223847;-0.431165761479662;-0.838135323053973;-0.260452777885674;0.868391080568684;0.106574186117463;-0.558464146596938;0.813736525946630;0.694591117684099];

Results.LCMXL = LCMXL(INPUT,Results,EstimOpt,OptimOpt);

xlswrite('LCMXL_results_demo.xls', Results.LCMXL.R_out)


%% ****************************     MIMIC     ****************************


% EstimOpt.StrNorm = ones(1,9); % standarize mean and variance of explanatory variables in INPUT.Xstr (default = true)

% Results.MIMIC.b0 = [0.148804376036183;-0.0364082689222229;-0.011064187742934;-0.0328575964505152;0.0424262689999236;0.00273733107195817;0.000235868427034212;-0.0768321543374594;-0.0720384450063044;-0.182575092847832;-0.0889680242634086;0.0379016122905587;-0.0244595986981180;-0.102660615070606;0.0539610288365757;0.0199342038767283;0.0718078258319252;0.0584470737842732;-0.201257402822855;0.0764346466677323;-0.163835971369101;-0.459578451640330;-0.547742335961172;-0.0378517363338953;0.131870616625717;0.157282552503199;0.223384083050719;-0.646448862649465;1.76638999126690;0.170229689466746;-0.467601747733821;-0.385812150452167;0.0639992488021068;0.354924420649614;-0.884183959319061;-0.456373780888765;-1.72488330525692;-0.819176767599829;-3.53925172461022;-0.614934148806317;-1.02462552933921;-0.124463517158109;-0.817505342443083;0.353472442572121;-1.44476685199243;0.294421923674619;-0.783669894549255;0.422497373690311;-2.23158336924593;0.576352808306931;-0.615136705928837;0.288301309162005;5.64281169897769;-2.24959315348704;0.292532147431253;-5.80848898507390;-2.28905990132861;-1.51817407957309;-1.95292496115563;0.668875562745533;-1.14733027425654;-0.158305917441182;0.371733658905934;-0.00971089741677728];
% Results.MIMIC = MIMIC(INPUT,Results,EstimOpt,OptimOpt);

% xlswrite('MIMIC.xls', Results.MIMIC.R_out)


%% ****************************     HMNL     ****************************


% EstimOpt.NLatent = 1; % number of Latent Variables
% EstimOpt.MeaMatrix = [] % Provide a matrix associating measurement equations with LV 
% EstimOpt.MeaExpMatrix = [] % Provide a matrix associating additional explanatory variables with measurement equations
% EstimOpt.MeaSpecMatrix = [2 2 1 2 2 0 5 2]; % Specify measurement equation type - 0 = OLS, 1 = MNL, 2 = ordered probit, 3 = Poisson, 4 = Negative Binomial, 5 = ZIP, 6 = ZINB
% EstimOpt.CountCut = 80 % Provide limit for censoring count variables (in XMea)
% EstimOpt.ScaleLV = 1; % Scale interacted with LV

Results.HMNL.b0 = [-1.06185535043816;0.347296925261114;0.535861862307341;0.00806290962208581;-0.0354294588168786;0.966640180671367;2.31437323890006;-0.416307402041866;-0.438947672859781;-0.527608417502472;-0.307608094813139;0.360015185154331];
 
Results.HMNL = HMNL(INPUT,Results,EstimOpt,OptimOpt);
  

%% ****************************     HLC     ****************************


Results.HLC.b0 = [1.19312701565716;0.281573204925181;0.537456981905208;-0.305421344236443;-0.183138203475598;1.10016741127977;-1.41314254530083;0.557004796772487;0.735144922319455;0.327296058102634;0.147110857939958;0.688667720913533;-1.46617998445493;5.18881748972799];

Results.HLC = HLC(INPUT,Results,EstimOpt,OptimOpt);

xlswrite('HLC_results_demo.xls', Results.HLC.R_out)


%% ****************************     HMXL     ****************************


EstimOpt.FullCov = 0;
 
Results.HMXL_d.b0 = [-1.06385396719045;0.502225968137010;0.762919876291738;0.0754979228406100;-0.00525178571432824;1.21463207189101;0.827435619057455;0.858379313513640;0.850626040781443;0.544743011345145;-0.231951625486396;1.04919002940943;2.76562009995691;-0.423410460175780;-0.475420579953603;-0.655841777057070;-0.454305416995704;-0.173435233502941];

Results.HMXL_d = HMXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     HMXL     ****************************
 

EstimOpt.FullCov = 1;
 
Results.HMXL.b0 = [-1.06385396719045;0.502225968137010;0.762919876291738;0.0754979228406100;-0.00525178571432824;1.21463207189101;0.827435619057455;0.858379313513640;0.850626040781443;0.544743011345145;-0.231951625486396;1.04919002940943;2.76562009995691;-0.423410460175780;-0.475420579953603;-0.655841777057070;-0.454305416995704;-0.173435233502941];

Results.HMXL = HMXL(INPUT,Results,EstimOpt,OptimOpt);

