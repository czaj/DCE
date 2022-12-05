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

DATA.Xm = [DATA.AGE==3, DATA.AGE==4, DATA.AGE==5]; 
EstimOpt.NamesM = {'e.g. 26-44','e.g. 45-64','e.g. 65+'}; 

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
INPUT.Xm = DATA.Xm(DATA.filter,:); % Explanatory variables of random parameters means / inreactions
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

B_backup = [-0.0649871757079609;0.516086083936041;0.764495593032389;0.175803035783724;0.0423877876196776;0.656304581992665;0.328582977085665;0.0112058878541658;-0.183656735546277;-0.0571270745882275;4.33728554883648e-05;-0.168997051214565;0.14460696302212;-0.0611137990713113;-0.204712689328378;0.169826650333787;0.10431372306787;0.0474591006391481;0.401634785987956;-0.0535136809608674;-0.218576712727205;-0.0601489890120321;0.00391361371887034;-0.20931481910203]; % provide starting values for the MNL model

Results.MNL = MNL(INPUT,Results,EstimOpt,OptimOpt);

B_backup = [0.157151789122682;0.486591518506085;0.608999367821945;0.187839637442765;0.0683659481870268;0.492009465386079;0;0;0;0;0;0.0631226149552788;0;0;0;0;0;0.153313208191371;0;0;0;0;0;0.0909225537519497]; % provide starting values for the MNL model

EstimOpt.BActive = [1;1;1;1;1;1;0;0;0;0;0;1;0;0;0;0;0;1;0;0;0;0;0;1];

Results.MNL = MNL(INPUT,Results,EstimOpt,OptimOpt);

%% ****************************     MXL_d     ****************************

% EstimOpt.FullCov = 0; % Specify if random parametrs are correlated or not (0 = default)
% EstimOpt.Dist = zeros(size(INPUT.Xa,2),1); % choose Xa parameters' distributions (0 - normal, -1 - constant, 1 - lognormal, 2 - Spike, 3 - triangular, 4 - Weibull, 5 - sinh-arcsinh, 6 - Johnson Sb, 7 - Johnson Su); 
% EstimOpt.Scores = 1; % Predict individual-specific paramter / WTP scores using bayes rule

EstimOpt.Dist = [0;0;0;0;0;1];

EstimOpt.BActive = ones(30);

B_backup = [-1.23685884591626;0.679719111089912;1.06706905559448;0.233782437696263;0.0670126109050347;-0.0947684440480008;3.28854699236889;0.86847941847826;0.88563006752372;0.631869482261344;-0.265050188044538;1.05668467646537;0.598254738657584;0.0224204349313097;-0.19815409257359;-0.0540320181435001;0.00397739850289974;-0.136819869420428;-0.0537098579662287;-0.0725358644643713;-0.273618863194524;0.244079388341653;0.165361901927099;0.243590109184194;0.600681298692093;-0.0938455827541106;-0.288702111002723;-0.0584351324695277;0.00873834035928742;-0.287535655048998]; % provide starting values for the MXL_d model

Results.MXL_d = MXL(INPUT,Results,EstimOpt,OptimOpt);

EstimOpt.BActive = [1;1;1;1;1;1;1;1;1;1;1;1;0;0;0;0;0;1;0;0;0;0;0;1;0;0;0;0;0;1];

B_backup = [-0.998340558566231;0.641933580167909;0.876322347471332;0.263578807426884;0.113723659256506;-0.168068608874872;3.27151058592506;0.866636798313879;0.891813015024776;0.636214264088218;-0.189404720947657;1.08566583709384;0;0;0;0;0;-0.0570679813701219;0;0;0;0;0;0.23687414990181;0;0;0;0;0;-0.0771831736781015]; % provide starting values for the MXL_d model

Results.MXL_d = MXL(INPUT,Results,EstimOpt,OptimOpt);

%% ****************************     MXL     ****************************

EstimOpt.FullCov = 1;

EstimOpt.BActive = ones(45);

B_backup = [-1.49322131885424;0.69563236031811;1.10135103764226;0.123264944924421;-0.026277952625769;-0.0865755832090546;3.21087225850338;-0.596202979217673;-0.49063429034651;-0.47730000852456;-0.35338219606974;-0.499318214435402;0.920806476561121;0.478462615773341;-0.0231308372989304;-0.227263684520831;0.272432645701717;0.907445407354726;0.304272354876634;0.145269456525003;-0.0926077485862265;0.75978570572526;0.418589139825705;-0.376890367099038;0.0137455342006848;-0.408611415627389;0.833128665825497;0.840920965533848;-0.0597812950922486;-0.324429699282103;-0.212297518209787;-0.067148729324853;-0.350311248021963;0.0320580978077737;-0.182554357504307;-0.392966646104291;0.117475935907562;0.096553135707362;0.0579761274143145;0.70715318965839;-0.261215078680759;-0.482579934676955;-0.221876702468782;-0.0711271489301003;-0.578831577339254]; % provide starting values for the MXL model

Results.MXL = MXL(INPUT,Results,EstimOpt,OptimOpt);

EstimOpt.BActive = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;0;0;0;0;0;1;0;0;0;0;0;1;0;0;0;0;0;1];

B_backup = [-1.29718514500468;0.607424270054228;0.821008037354993;0.0545687717354477;-0.0562189140430507;-0.301590524006771;3.09756129099751;-0.637320064110808;-0.477804890994866;-0.569882494138378;-0.401182103926001;-0.383663359623434;0.827191030269984;0.315500190605376;-0.312278351296585;-0.414230379742765;0.583879532191886;0.873648988119204;0.20221762454249;0.0595021667608166;-0.0082894647960951;0.538933209248466;0.0988976069338563;-0.0951799304213314;-0.00727210872870863;0.593291746219842;0.809919518324464;0;0;0;0;0;-0.110488481297237;0;0;0;0;0;0.383308086962723;0;0;0;0;0;-0.0629454473135256]; % provide starting values for the MXL model

Results.MXL = MXL(INPUT,Results,EstimOpt,OptimOpt);
