% DCE demo script - 2022-12-05
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


EstimOpt.DataFile = ('NEWFOREX_DCE_demo1.mat'); % specify the name of the data file
DATA = load(EstimOpt.DataFile); % load data, store it in a structure DATA


%% ****************************  data description  ****************************

tableDCE = [DATA.ID, DATA.ALT, DATA.CT, DATA.Choice, DATA.SQ, DATA.GOS, DATA.CEN, DATA.VIS==2, DATA.VIS==1, -DATA.FEE/4/10, DATA.TIME, DATA.AGE];
tableDCE = array2table(tableDCE, 'VariableNames', {'Id', 'Alt', 'CS', 'Choice', 'Status quo', 'Passive protection of all commertial forests', 'Passive protection of all natural regeneration forests', 'Reducing the number of visitors to 7500 per day', 'Reducing the number of visitors to 5000 per day', '-Cost (10 EUR)', 'Time', 'Age'});

%% ****************************  model specification  ****************************

EstimOpt.ProjectName = 'demo'; % name of the project used for naming xls results files

DATA.Y = DATA.Choice; % exlpanained variable (choice indicator)

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

% Optional: pass starting values either using B_backup (global variable automatically saved at each iteration) or Results.MNL.b0
B_backup = [0.15995426181761282503;0.48658845316475168863;0.60908549594851302267;0.18827611825172022031;0.068210544286035651451;0.56989320442201241157];

Results.MNL = MNL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     MXL_d     ****************************

% EstimOpt.FullCov = 0; % Specify if random parametrs are correlated or not (0 = default)
% EstimOpt.Dist = zeros(size(INPUT.Xa,2),1); % choose Xa parameters' distributions (0 - normal, -1 - constant, 1 - lognormal, 2 - Spike, 3 - triangular, 4 - Weibull, 5 - sinh-arcsinh, 6 - Johnson Sb, 7 - Johnson Su); 
% EstimOpt.Scores = 1; % Predict individual-specific paramter / WTP scores using bayes rule

EstimOpt.Dist = [0;0;0;0;0;1];

B_backup = [-1.0862014791638248123;0.64232140635357226621;0.86507031518939914161;0.26073942118058013317;0.10734913592046885222;-0.099211061936221894841;3.0447049702322388498;0.87224682281311016752;0.89381609008317963738;0.61181029893998228886;-0.23653153447499053463;1.2642393272248089175];

Results.MXL_d = MXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     MXL     ****************************

EstimOpt.FullCov = 1;

B_backup = [-1.1972997845872539457;0.56466212013198713304;0.84065147282087948621;0.056448091884404236196;-0.042022520905620995568;-0.24953593149086314429;3.1494409625688204457;-0.55420400070170305895;-0.46052293523460213764;-0.52370556394074641027;-0.38867189226771919897;-0.44406792688772656064;0.89914624166189760501;0.43150921427762245486;-0.097950311613189844362;-0.24948995878247606783;0.19237824462260666447;0.90826690774189344779;0.23367547499733645755;0.09495176021286447221;0.047413224981750615172;0.72531243869184869322;0.33234481455619935275;-0.3466234016364158621;-0.23534646202378431412;0.33581892952799974328;0.89817829797193338148];

Results.MXL = MXL(INPUT,Results,EstimOpt,OptimOpt);

%% ****************************     GMXL_d     ****************************


EstimOpt.FullCov = 0;
% EstimOpt.Gamma0 = 0.5; % 0 and 1 impose GMXL type II and I, respectively. Value in between allows for gamma to be freely estimated.

Results.GMXL_d.b0 = [-1.7108863135595624438;0.98831747151039406329;1.2935707690084969901;0.46209169702299801585;0.20968570810056438858;0.37793010995203574209;4.4673095194292340437;1.043229948524629469;1.0984008868841879103;0.70670541206687353952;-0.39023578534185354716;1.2519126775641125082;1.0231506638682434929;0];

Results.GMXL_d = GMXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     GMXL     ****************************


EstimOpt.FullCov = 1;

Results.GMXL.b0 = [-1.7852911584731743222;0.72373528558407485001;1.2377118874761519063;-0.099115785911607579006;-0.12584093376937130482;0.22450807200333830482;4.1567107205111328838;-0.70096944363116819865;-0.35766720154211589788;-0.78182573115662601371;-0.63022664945990447549;-0.29603513431704303605;1.0389857154922588212;0.32968282866545717269;-0.33208095602767295773;-0.35785032597439314639;-0.034123762610221146374;1.0425516581650573489;0.30963134914172302237;0.12826542670548937708;-0.039161658056455932175;0.75223048547544346665;0.31619673539266845985;-1.0606403878294381471;0.23322962287703888351;-0.21307644297545774714;0.010956608732913564533;0.9457411504244138678;0];

Results.GMXL = GMXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     LC     ****************************


% EstimOpt.NClass = 2; % number of latent classes
% EstimOpt.ClassScores = 1; % Predict individual-specific class membership probabilities using Bayes rule
% EstimOpt.BActiveClass = ones(EstimOpt.NVarA,1); % use 0 to impose coefficient equality for all latent classes

Results.LC.b0 = [1.1931370757907797664;0.28157647154191456362;0.53747172960122202579;-0.30543839419643020738;-0.18314773628703553965;1.1001689372825087521;-1.4131343817410375596;0.55700292525106476216;0.73514118772222236675;0.32729495687496545919;0.14711082672045694419;0.68866433413028960153;-0.42772878402315994695];

Results.LC = LC(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     LCMXL_d     ****************************


EstimOpt.FullCov = 0;
% EstimOpt.Scores = 1; % Predict individual-specific paramters / class membership probabilities using Bayes rule

Results.LCMXL_d.b0 = [0.038080742683099931545;0.22818487322238015236;0.31742761286104159701;-0.12565440286973314499;-0.16235555344963953361;-0.32269702082043832947;-3.8239003331275562836;1.4181961167717460626;1.9364405524738381725;0.97610731083980306622;0.60158701069160391839;0.39257775070026945663;2.0568840041824745235;0.41568420963400160018;-0.32521857358722161546;0.39496011085925880613;0.33549619261823931948;1.3919258768556674877;2.1335591515587712941;1.2754360293186153275;1.241915301363331281;-0.86532199582279778483;-0.12129102879446600205;1.0419217222020118463;0.4290906748324989084];

Results.LCMXL_d = LCMXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     LCMXL     ****************************


EstimOpt.FullCov = 1;

Results.LCMXL.b0 = [-0.49533744048379146907;0.3345535558943508736;0.43933760416634476398;-0.025983459552296517964;-0.098375424360727506401;-0.50914554382368215624;-4.6243064525101669204;1.8501737079220275106;2.4214711418524972331;0.26017690840050733403;0.22445115381500660434;0.91191758572236814029;2.3748920898452947625;-0.25789048366503941612;-0.096815829957440494025;-0.18922068541960601618;-0.11731695569289021797;-0.36768176821739206872;0.46555361664951089296;0.46004749263418431848;-0.0089295628885870444863;-0.31736858792622141268;0.61589581626025191596;-0.088192768628825149446;-0.50680597548561934218;-0.20653122192894951548;0.58648495413701107193;0.050389513562356540166;-0.076054432045359915415;0.065979284543692146014;-0.063022683357202707866;0.078233805234784084548;-0.65115625175552072079;3.164436905354462759;-0.550106186922620366;0.44621839092480886091;-0.96300575020698830908;-0.74218949398927658301;-0.48395585880685598745;1.843701090019083999;0.57992594421850174324;-0.406089529642018221;-0.53084048603643796405;-0.077371370189640589765;1.5859108264662442611;0.80759194139771495191;0.78729258063985496641;-0.45736922211237379665;-1.3058354023013560852;-0.53963199128310179731;0.95746007689095613546;-0.084244581825056383262;-0.67882402540494679588;1.2516136203775070079;0.44015894372226643805];

Results.LCMXL = LCMXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     MIMIC     ****************************


% EstimOpt.StrNorm = ones(1,9); % standarize mean and variance of explanatory variables in INPUT.Xstr (default = true)

Results.MIMIC.b0 = [2.259614998884304429e-07;0.026333461635994544203;-0.00036443559706179472907];

Results.MIMIC = MIMIC(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     HMNL     ****************************


% EstimOpt.NLatent = 1; % number of Latent Variables
% EstimOpt.MeaMatrix = [] % Provide a matrix associating measurement equations with LV 
% EstimOpt.MeaExpMatrix = [] % Provide a matrix associating additional explanatory variables with measurement equations
% EstimOpt.MeaSpecMatrix = [2 2 1 2 2 0 5 2]; % Specify measurement equation type - 0 = OLS, 1 = MNL, 2 = ordered probit, 3 = Poisson, 4 = Negative Binomial, 5 = ZIP, 6 = ZINB
% EstimOpt.CountCut = 80 % Provide limit for censoring count variables (in XMea)
% EstimOpt.ScaleLV = 1; % Scale interacted with LV

Results.HMNL.b0 = [-1.0642222680313289107;0.34377512643953744842;0.52930541202094094633;0.0071248707849624007965;-0.037349620942344494146;-0.063724804011940666681;-2.1967769270033246087;0.42363916941303325636;0.45528441671209368691;0.52874334552423507549;0.31621152326530865828;-0.44236918577681716425;-0.00012538777785397854582;0.1701197868570860916;-0.01470207763384337625];
 
Results.HMNL = HMNL(INPUT,Results,EstimOpt,OptimOpt);
  

%% ****************************     HLC     ****************************


Results.HLC.b0 = [-1.4149022250171274795;0.5572779698668448578;0.73555131099661508198;0.32697389325297043783;0.14690540003923541201;0.68897251915770352237;1.1951283179144924773;0.28069271706284110168;0.53438472226581745783;-0.30100045172269401794;-0.18109666466492840842;1.0962545518833841651;0.52197311822648784219;1.032125214689078696;6.5613580294080172989e-05;0.32228133084264731778;-0.054813545186653003793];

Results.HLC = HLC(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     HMXL     ****************************


EstimOpt.FullCov = 0;
 
Results.HMXL_d.b0 = [-1.2435718117557386098;0.62331054872172608761;0.85937724623059530416;0.2648166947029513274;0.12226274957204133487;-0.1851837208291867154;2.6578400043858896318;0.7629684931701193884;0.66814232268017037519;0.60194898849233657856;-0.33757846049636686114;1.1626314688056507141;-1.656050506481095308;0.62888540215818766743;0.8260252549495805674;0.32030504576832879993;0.1342507115915810012;0.52377759905227094794;0.001149170410149287936;0.38103374084134006283;-0.078703987343929790454];

Results.HMXL_d = HMXL(INPUT,Results,EstimOpt,OptimOpt);


%% ****************************     HMXL     ****************************
 

EstimOpt.FullCov = 1;
 
Results.HMXL.b0 = [-1.4133841505896469481;0.57652306496438288086;0.88815958411205353507;0.0010762528533418815133;-0.07655466699069995995;-0.26244840358225901555;2.501571925443574429;-0.30219795626342754735;0.022048458946810175346;-0.44225878102294646776;-0.375985992631937882;-0.42719477629250796058;0.84354052001164570029;0.21116062525503651037;-0.1052554962846255493;-0.23358149070720762852;-0.039381527158466607397;0.76966494081380332215;0.29727126190010716433;0.14500279414232605801;-0.012844590982789004172;0.75263935453897723438;0.34403878959320588482;-0.78561915753581479382;0.21788126951152950173;-0.86797094578113431762;0.25377798482872326868;-1.8904615232394161861;0.59073331663885375598;0.76932003816915550676;0.26228449857165836434;0.10312563413344265351;0.50390085112316074234;-0.00077990497001570197738;0.44766547728213551549;-0.11135022527092809652];

Results.HMXL = HMXL(INPUT,Results,EstimOpt,OptimOpt);