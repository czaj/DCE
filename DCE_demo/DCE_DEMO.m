
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

DATA.TIME = zeros(size(DATA.TIME_1));
for i = 1:36
    DATA.TIME = DATA.TIME + eval(['DATA.TIME_',num2str(i)]);
end
DATA.TIME = VarNorm(DATA.TIME);

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

DATA.Xmea = [DATA.TIME]; % 
EstimOpt.NamesMea = {'DCE time (normalized)'};

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
INPUT.Xmea = DATA.Xmea(DATA.filter,:); % Measurment equations variables (Hybrid models)
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
 
Results.MXL_d.b0 = [-0.729340119066748;0.633825526889206;0.861187549808267;0.267036729140991;0.121825890107902;1.21950771816509;3.27178757659456;0.872309731140975;0.905806187595547;0.588931734234715;-0.233694908995634;1.06279280385645]; % provide starting values for the MXL_d model

Results.MXL_d = MXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     MXL     ****************************


EstimOpt.FullCov = 1;

Results.MXL.b0 = [-1.07526723050544;0.543123503778332;0.834308099998921;0.0324671780476724;-0.0272456319736574;1.22152638053686;2.95344147563694;-0.494625407376323;-0.422197186338681;-0.542075443286223;-0.367172207868796;-0.274374550958867;0.896088845617134;0.377570410797849;-0.110142017955691;-0.211424687654118;0.0150290537518090;0.941065837069406;0.235676071491103;0.115643919845501;0.0363027067083060;0.705842268512163;0.363638702608856;-0.490593045318898;-0.327832993803731;0.208727891050920;0.909216688885196]; % provide starting values for the MXL model

Results.MXL = MXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     GMXL_d     ****************************


EstimOpt.FullCov = 0;
% EstimOpt.Gamma0 = 0.5; % 0 and 1 impose GMXL type II and I, respectively. Value in between allows for gamma to be freely estimated.

Results.GMXL_d.b0 = [-1.19814256548968;0.937481653339174;1.28041595724521;0.454247168674383;0.224605301498210;1.83979189784192;4.63801659339025;1.05870684175096;1.02438521536684;0.658168936363674;-0.294334044502617;1.59695735632123;0.951490333930232;0]; % provide starting values for the GMXL_d model

Results.GMXL_d = GMXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     GMXL     ****************************


EstimOpt.FullCov = 1;

Results.GMXL.b0 = [-1.49564358686303;0.754208806660050;1.24408224551453;-0.0322270839069335;-0.0724990493564693;1.79215431042702;3.98383226312750;-0.548938365209824;-0.326144589187290;-0.805109692589865;-0.586632842562863;0.0534460380173224;1.05752902315142;0.298097162930150;-0.288543223930740;-0.294796627078254;-0.119627530598077;1.00859947894374;0.123787344243268;0.0469044121704953;0.0426217195246578;0.558485565843005;0.266512052917830;-1.15013213284956;0.0105333844410333;-0.0117191362785322;0.556100451454097;0.932055063730309;0]; % provide starting values for the GMXL model

Results.GMXL = GMXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     LC     ****************************


% EstimOpt.NClass = 2; % number of latent classes
% EstimOpt.ClassScores = 1; % Predict individual-specific class membership probabilities using Bayes rule
% EstimOpt.BActiveClass = ones(EstimOpt.NVarA,1); % use 0 to impose coefficient equality for all latent classes

Results.LC.b0 = [1.19313730540666;0.281576712159116;0.537472373318735;-0.305438503065379;-0.183148111515675;1.10016928328555;-1.41313375267448;0.557003045472699;0.735140968672557;0.327295416756785;0.147111072049196;0.688664095386527;-0.427725860586280];

Results.LC = LC(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     LCMXL_d     ****************************


EstimOpt.FullCov = 0;
% EstimOpt.Scores = 1; % Predict individual-specific paramters / class membership probabilities using Bayes rule

Results.LCMXL_d.b0 = [0.409558531584716;0.270570487898003;0.339643271373607;-0.115056580801768;-0.169857260397970;1.03239449212084;-3.65424085051403;1.38438510740756;1.97737573894267;1.04246944408380;0.654559710856222;1.71859421544641;2.34650344954097;0.437253695346915;-0.280352959440886;0.323874913988266;0.370359174059547;0.872235035826814;1.80348362867858;1.32238896275705;1.28458067725295;-0.869833667411202;0.0723471798054098;1.31214790027106;0.549239831186432];

Results.LCMXL_d = LCMXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     LCMXL     ****************************


EstimOpt.FullCov = 1;

Results.LCMXL.b0 = [0.346789921894470;0.248012160717817;0.433353114857925;-0.192321423511116;-0.210562175403753;0.921407703284509;-3.99828820109820;1.53700973228666;2.39595388630782;0.649489515501281;0.507850982258279;2.17297313750519;2.33844322054203;-0.184585632885740;0.115896416513661;-0.257912979969864;-0.110999674529241;0.150980413176492;0.491769290653490;0.421416138495307;0.176924628688469;-0.195771091850034;0.222285865046263;0.0869510739916477;-0.665501852284466;-0.504082499919503;0.473251619733652;-0.0298299927446039;-0.235738181758925;0.128888149519537;0.0560011863835082;0.338275005344323;-0.391700195954950;2.17674399856972;-0.0675018624774421;0.611221657655554;-0.619409681663272;-0.351587322295996;-0.476248026311962;1.55735397261548;0.240122901277161;-0.303275494327079;-0.375909481351684;-0.327549712832126;1.10650543967390;0.666082786615312;0.376834171061103;-0.396376297675030;-0.796850074559161;-0.184985629097956;0.793422239243634;0.188264822815658;-0.547350816192404;1.05671002215609;0.638065578657326];

Results.LCMXL = LCMXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     MIMIC     ****************************


% EstimOpt.StrNorm = ones(1,9); % standarize mean and variance of explanatory variables in INPUT.Xstr (default = true)

Results.MIMIC.b0 = [2.40087403627601e-07;0.0263334597565126;-0.000364490204859268];

Results.MIMIC = MIMIC(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     HMNL     ****************************


% EstimOpt.NLatent = 1; % number of Latent Variables
% EstimOpt.MeaMatrix = [] % Provide a matrix associating measurement equations with LV 
% EstimOpt.MeaExpMatrix = [] % Provide a matrix associating additional explanatory variables with measurement equations
% EstimOpt.MeaSpecMatrix = [2 2 1 2 2 0 5 2]; % Specify measurement equation type - 0 = OLS, 1 = MNL, 2 = ordered probit, 3 = Poisson, 4 = Negative Binomial, 5 = ZIP, 6 = ZINB
% EstimOpt.CountCut = 80 % Provide limit for censoring count variables (in XMea)
% EstimOpt.ScaleLV = 1; % Scale interacted with LV

Results.HMNL.b0 = [-1.06196392092922;0.347344784715878;0.532317757337924;0.0104716796358509;-0.0350388762293367;0.963524522005476;-2.32452122609725;0.418410463990034;0.448621020658069;0.525621052323958;0.309153772746186;-0.349884134351600;-7.54965570514627e-05;0.168937364373293;-0.0145012986224121];
 
Results.HMNL = HMNL(INPUT,Results,EstimOpt,OptimOpt);
  

%% ****************************     HLC     ****************************


Results.HLC.b0 = [-1.41490297455231;0.557277757762933;0.735551523100527;0.326974529564706;0.146904551623588;0.688973155469439;1.19512814044071;0.280692929166753;0.534383873850170;-0.301000451722694;-0.181095816249281;1.09625319751800;0.521972906122576;1.03212542019115;6.49772685585255e-05;0.322281755050471;-0.0548131209788293];

Results.HLC = HLC(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     HMXL     ****************************


EstimOpt.FullCov = 0;
 
Results.HMXL_d.b0 = [-1.06385396719045;0.502225968137010;0.762919876291738;0.0754979228406100;-0.00525178571432824;1.21463207189101;0.827435619057455;0.858379313513640;0.850626040781443;0.544743011345145;-0.231951625486396;1.04919002940943;2.76562009995691;-0.423410460175780;-0.475420579953603;-0.655841777057070;-0.454305416995704;-0.173435233502941];

Results.HMXL_d = HMXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     HMXL     ****************************
 

EstimOpt.FullCov = 1;
 
Results.HMXL.b0 = [-1.06385396719045;0.502225968137010;0.762919876291738;0.0754979228406100;-0.00525178571432824;1.21463207189101;0.827435619057455;0.858379313513640;0.850626040781443;0.544743011345145;-0.231951625486396;1.04919002940943;2.76562009995691;-0.423410460175780;-0.475420579953603;-0.655841777057070;-0.454305416995704;-0.173435233502941];

Results.HMXL = HMXL(INPUT,Results,EstimOpt,OptimOpt);

