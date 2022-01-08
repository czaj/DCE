%% Demo #1 for MDCEV
        % Basic simulation for Alpha and Gamma profiles
        % B denotes assumed coefficients for attributes (Xa)
        % Alpha and Gamma are coefficients for each profile
        % Scale is a std. dev. of error term
clear
clc
rng(10000001) % Setting a random seed
global B_backup;

EstimOpt.NCT = 1;
EstimOpt.NP = 1000;
EstimOpt.NAlt = 5;

OptimOpt = optimoptions('fminunc');
%OptimOpt.Algorithm = 'quasi-newton';
OptimOpt.GradObj = 'on';
% OptimOpt.Hessian = 'user-supplied'; % ('off'), only used by trust-region
OptimOpt.MaxIter = 1000;
OptimOpt.FunValCheck = 'on';
OptimOpt.Diagnostics = 'on';
OptimOpt.MaxFunEvals = 1e6; %Maximum number of function evaluations allowed (1000)
OptimOpt.OutputFcn = @outputf;
OptimOpt.TolFun = 1.e-7;
OptimOpt.TolX = 1.e-7;

N = EstimOpt.NP*EstimOpt.NCT;

EstimOpt.Profile = 2; % = 1 if Alpha-profile, = 2 is Gamma

%% Generating data
Price = (1+poissrnd(8, [EstimOpt.NAlt, N]))/4; % Price vec
Z = unifrnd(0,1, [EstimOpt.NAlt, N]); % Quality
ASC = [ones(1,N); zeros(4,N)]; % ASC for first alternative
Income = 1+ lognrnd(2, 0.6, [1000,1]);
eps = unifrnd(0,1, [EstimOpt.NAlt, N]); % error term
eps = -log(-log(eps)); % EV distribution

B = [-1; 1]; 
Scale = 1;
Phi = exp(B(1)*ASC + B(2)*Z + Scale*eps); 
if EstimOpt.Profile == 1
    Alpha = unifrnd(-1.5, -0.1, [EstimOpt.NAlt, 1]);
    Gamma = ones(EstimOpt.NAlt, 1); % Gamma coefs
else
    Alpha = zeros(EstimOpt.NAlt, 1); %  
    Gamma = unifrnd(0.5, 1.5, [EstimOpt.NAlt, 1]);
end

%% Generating dependent variable

Y = NaN(EstimOpt.NAlt, N);
for i = 1:N
    Start = Income(i)/EstimOpt.NAlt;
    Start = Start./Price(2:end,i); % Starting values for optimization
   funx = @(x) -utilityMDCEV(x,Price(:,i), Phi(:,i), Alpha, Gamma, Income(i), EstimOpt.Profile);
   % Restrictions: greater than 0 and less than income
   A = -eye(EstimOpt.NAlt-1);
   A = [A; Price(2:end,i)'];
   b = [zeros(EstimOpt.NAlt-1,1); Income(i)];
   tmp = fmincon(funx,Start,A,b);
 %   tmp = fminunc(funx, log(Start));
   % x = exp(tmp);
    x = tmp;
    x1 = (Income(i) - sum(Price(2:end,i).*x,1))/Price(1,i);
    Y(:,i) = [x1; x];
end

%% Data transformation

INPUT.Xa = [reshape(ASC, [EstimOpt.NAlt*N,1]), reshape(Z, [EstimOpt.NAlt*N,1])];
INPUT.Y = reshape(Y, [EstimOpt.NAlt*N,1]);
INPUT.I = Income;
INPUT.priceMat = reshape(Price, [EstimOpt.NAlt*N,1]);
INPUT.Y(INPUT.Y < 0.0001) = 0; % Needed?

% Test: how many goods individuals choose
% A = Y;
% AA = sum(A > 0.0001)';
% [mean(AA == 1), mean(AA == 2), mean(AA == 3), mean(AA == 4), mean(AA == 5)] 
% pause;

%% Estimation
if EstimOpt.Profile == 1
    B0 = [B; Alpha; Scale]; % True values of parameters (need to be transformed for starting values)
else
    B0 = [B; Gamma; Scale]; % True values of parameters
end

B_backup = [];
EstimOpt.HessEstFix = 1;
Results.MDCEV = MDCEV(INPUT,EstimOpt,OptimOpt);




%% Calculate utility (additional function)
function U = utilityMDCEV(x,Price, Phi, Alpha, Gamma, I, Prof)
   x1 = (I - sum(x.*Price(2:end),1))/Price(1);
   x = [x1; x];
   if Prof == 1 % Alpha
       U = (x./Gamma + 1).^Alpha - 1;
       U = Gamma.*U./Alpha;
   else % Gamma
       U = log(x./Gamma + 1);
       U = Gamma.*U;
   end
   U = sum(Phi.*U,1);
end