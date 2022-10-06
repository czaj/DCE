function [D, Err] = MDCEV_demand(Results, num)

% Functions takes as an argument Results from MDCEV estimation
    % num is empty if there is no numeraire
    % otherwise num should indicate which alternative is numeraire

% D is an NAlt x NCT*NP matrix with predictions of Demand
%   Err is prediction error

% Example:
% [D, Err] = MDCEV_demand(Results.MDCEV, []);

errorx = 0.000001; % acceptable error for bisection (only alpha-profile)

%% Variables and parameters
EstimOpt = Results.EstimOpt;
if isfield(EstimOpt,'NSim') == 0
    EstimOpt.NSim = 500;
end
Xa = Results.INPUT.Xa;
Y = Results.INPUT.Y;
I = Results.INPUT.I'; % I should 1 x N
Xm = Results.INPUT.Xm; % covariates (socio-demographic) (NxNVarM)
priceMat = Results.INPUT.priceMat; % matrix of prices of each alternative (NAltxN)

variables = Results.bhat;
NVarA = EstimOpt.NVarA; % Number of attributes
NVarP = EstimOpt.NVarP; % Number of parameters for alpha/gamma profile
NVarM = EstimOpt.NVarM;
Profile = EstimOpt.Profile; % Utility function version
SpecProfile = EstimOpt.SpecProfile;
NAlt = EstimOpt.NAlt; % Number of alternatives
N = size(Y,2); % Number of decisions

% define variables to be optimized
betas = variables(1:NVarA); % betas
b_profile = variables(NVarA+1:NVarA+NVarP*(1+NVarM));
scale = exp(variables(NVarA+NVarP*(1+NVarM)+1)); % scale parameter (sigma); one for all dataset

%% Specifying alphas and gammas
if Profile == 1 % alpha profile
    if NVarM == 0
        alphas = b_profile;
        alphas = 1 - exp(-alphas); % between 0 and 1
    else
        alpha = reshape(b_profile, [NVarM+1, NAlt]);
        alphas = 1 - exp(-Xm*alpha)'; % NAlt x N
    end
    gammas = ones(size(alphas));
elseif Profile == 2 % gamma profile
    if NVarM == 0
        gammas = b_profile;
        gammas = exp(gammas); % greater than 0
    else
        gamma = reshape(b_profile, [NVarM+1, NAlt]);
        gammas= exp(Xm*gamma)'; % NAlt x N
    end
    alphas = zeros(size(gammas));
elseif Profile == 3 % SpecProfile
    indx = length(unique(SpecProfile(1,SpecProfile(1,:) ~= 0))); % unique alphas
    if NVarM == 0
        a = b_profile(1:indx);
        a = 1 - exp(-a);
        a = a(SpecProfile(1,SpecProfile(1,:) ~= 0)); 
        alphas = zeros(NAlt,1);
        alphas(SpecProfile(1,:) ~= 0) = a;
    
        g = b_profile(indx+1:end);
        g = exp(g); % greater than 0
        g = g(SpecProfile(2,SpecProfile(2,:) ~= 0)); 
        gammas = ones(NAlt,1);
        gammas(SpecProfile(2,:) ~= 0) = g;
    else

        a = reshape(b_profile(1:indx*(NVarM+1)), [NVarM+1, indx]);
        alpha = 1 - exp(-Xm*a)';
        alpha = alpha(SpecProfile(1,SpecProfile(1,:) ~= 0),:); 
        alphas = zeros(NAlt,N);
        alphas(SpecProfile(1,:) ~= 0,:) = alpha;
    
        g = reshape(b_profile(indx*(NVarM+1)+1:end), [NVarM+1, NVarP - indx]);
        gamma = exp(Xm*g)'; % greater than 0
        gamma = gamma(SpecProfile(2,SpecProfile(2,:) ~= 0),:); 
        gammas = ones(NAlt,N);
        gammas(SpecProfile(2,:) ~= 0,:) = gamma;
    end
end
 
%% Demand prediction

betasZ = reshape(Xa*betas(1:NVarA), [NAlt, N]); % beta*X part of utility
% generating error terms
hm1 = sobolset(NAlt,'Skip',1,'Leap',0);
hm1 = scramble(hm1,'MatousekAffineOwen');
eps = net(hm1,N*EstimOpt.NSim); % this takes every point:
clear hm1;
eps = reshape(-log(-log(eps')),[NAlt, EstimOpt.NSim, N]);
eps = permute(eps, [1 3 2]);
MU = exp(betasZ + scale*eps - log(priceMat)); % Price-adjusted baseline utility (NAlt x N x NSim)

pg = priceMat.*gammas;
inva = 1./(1-alphas);
 Vec = (1:NAlt)';
% The easier case - alphas constant across alternatives

if sum(SpecProfile(1,:)) == 0  || length(unique(SpecProfile(1,:)')) == 1
    MUa = MU.^inva;
    pgMUa = pg.*MUa;
    Dsim = zeros(NAlt, N, EstimOpt.NSim);
%     Lam = (I + sum(pg,1))./sum(pg.*MUa,1);
%     Dsim = (MUa.*Lam - 1).*gammas;
%     Lam = Lam.^(alphas-1);
%     NotChosen = Lam > MU;

%     sum(NotChosen(:))
    if isempty(num) % there is no numeraire
        for i = 1:N
            for j = 1:EstimOpt.NSim
                cond = 0;
                k = 2;
                % To start with the alternative with highest MU is consumed
                [MUsort,Indx] = sort(MU(:,i,j),'descend');
                Lam_k_tmp = (I(i) + pg(Indx(1),i))./pgMUa(Indx(1),i,j);
                Lam_k = Lam_k_tmp.^(alphas(Indx(1))-1); % This would need to be adjusted for individual specific alphas
                
                while cond == 0
                    if k > NAlt
                        cond = 1;
                    else
                        if Lam_k <= MUsort(k)
                            Lam_k_tmp = (I(i) + sum(pg(Indx(1:k),i),1))./sum(pgMUa(Indx(1:k),i,j),1);
                            Lam_k = Lam_k_tmp.^(alphas(Indx(1))-1);
                            k = k+1;
                        else
                            cond = 1;
                        end
                    end
                end
            % k-1 alternatyw kosumowanych
                Dsim(Indx(1:(k-1)),i,j) = (Lam_k_tmp.*MUa(Indx(1:(k-1)),i,j) - 1).*gammas(Indx(1:(k-1))); 
            end
        end
    else % numeraire
        if num ~= NAlt
            % put num at the end
            Vect = [Vec(Vec ~= num); num];
            pg = pg(Vect,:);
            pgMUa = pgMUa(Vect,:,:);
            alphas = alphas(Vect);
        end
        for i = 1:N
            for j = 1:EstimOpt.NSim

                cond = 0;
                k = 2;
                % To start with only numeraire is consumed
                Lam_k_tmp = (I(i) + pg(NAlt,i))./pgMUa(NAlt,i,j);
                Lam_k = Lam_k_tmp.^(alphas(NAlt)-1); % This would need to be adjusted for individual specific alphas
                [MUsort,Indx] = sort(MU(Vec ~= NAlt,i,j),'descend');
                
                while cond == 0
                    if k > NAlt
                        cond = 1;
                    else
                        if Lam_k <= MUsort(k-1)
                            FirstAlts = [Indx(1:(k-1)); NAlt];
                            Lam_k_tmp = (I(i) + sum(pg(FirstAlts,i),1))./sum(pgMUa(FirstAlts,i,j),1);
                            Lam_k = Lam_k_tmp.^(alphas(NAlt)-1);
                            k = k+1;
                        else
                            cond = 1;
                        end
                    end
                end
            % k-1 alternatyw kosumowanych
                Dsim(FirstAlts,i,j) = (Lam_k_tmp.*MUa(FirstAlts,i,j) - 1).*gammas(FirstAlts); 
            end
        end
    end
else % Alpha profile
   % error('Demand prediction works only if alphas are constant across utilities')
      MUa = MU.^inva;
    pgMUa = pg.*MUa;
    Dsim = zeros(NAlt, N, EstimOpt.NSim);

    if isempty(num) % there is no numeraire
        for i = 1:N
%             i
            for j = 1:EstimOpt.NSim
%                 j
                cond = 0;
                k = 1;
                % To start with only numeraire is consumed
                [MUsort,Indx] = sort(MU(:,i,j),'descend');
                Lam_hat = MU(Indx(2),i,j);
                E_hat = ehat(Lam_hat, MUa(Indx(1),i,j), inva(Indx(1)), pg(Indx(1),i));
                while cond == 0
                    if E_hat < I(i)
                        k = k+1; % 2 or more
                        if k < NAlt
                            Lam_hat = MU(Indx(k+1),i,j);
                            E_hat = ehat(Lam_hat, MUa(Indx(1:k),i,j), inva(Indx(1:k)), pg(Indx(1:k),i));
                        else
                            Lam_low = 0;
                            Lam_up = MU(Indx(k),i,j);
                            % go to bisection
                            lam = bisect(Lam_low, Lam_up, I(i), errorx, MUa(Indx(1:k),i,j), inva(Indx(1:k)), pg(Indx(1:k),i));
                            cond = 1;
                        end
                    else
                        Lam_low = MU(Indx(k),i,j);
                        Lam_up = MU(Indx(k+1),i,j);
                        lam = bisect(Lam_low, Lam_up, I(i), errorx, MUa(Indx(1:k),i,j), inva(Indx(1:k)), pg(Indx(1:k),i));
                        cond = 1;
                    end
                end
            % k-1 alternatyw kosumowanych
            Dsim(Indx(1:k),i,j) = (MUa(Indx(1:k),i,j).*(lam.^(-inva(Indx(1:k))))-1).*gammas(Indx(1:k));
%                 Dsim(FirstAlts,i,j) = (Lam_k_tmp.*MUa(FirstAlts,i,j) - 1).*gammas(FirstAlts); 
            end
        end
    else % numeraire
        if num ~= NAlt
            % put num at the end
            Vect = [Vec(Vec ~= num); num];
            pg = pg(Vect,:);
%             pgMUa = pgMUa(Vect,:,:);
%             alphas = alphas(Vect);
            inva = inva (Vect);
           
        end
        for i = 1:N
            for j = 1:EstimOpt.NSim

                cond = 0;
                k = 1;
                % To start with only numeraire is consumed
                [MUsort,Indx] = sort(MU(Vec ~= NAlt,i,j),'descend');
                Lam_hat = MU(Indx(1),i,j);
                E_hat = ehat(Lam_hat, MUa(NAlt,i,j), inva(NAlt), pg(NAlt,i));
                while cond == 0
                    if E_hat < I(i)
                        k = k+1; % 2 or more
                        FirstAlts = [Indx(1:(k-1)); NAlt];
                        if k < NAlt
                            Lam_hat = MU(Indx(k),i,j);
                            E_hat = ehat(Lam_hat, MUa(FirstAlts,i,j), inva(FirstAlts), pg(FirstAlts,i));
                        else
                            Lam_low = 0;
                            Lam_up = MU(Indx(k-1),i,j);
                            % go to bisection
                            lam = bisect(Lam_low, Lam_up, I(i), errorx, MUa(FirstAlts,i,j), inva(FirstAlts), pg(FirstAlts,i));
                            cond = 1;
                        end
                        k = k+1;
                    else
                        FirstAlts = [NAlt; Indx(1:(k-1))];
                        Lam_low = MU(FirstAlts(k),i,j);
                        Lam_up = MU(Indx(k),i,j);
                        FirstAlts = [Indx(1:(k-1)); NAlt];
                        lam = bisect(Lam_low, Lam_up, I(i), errorx, MUa(FirstAlts,i,j), inva(FirstAlts), pg(FirstAlts,i));
                        cond = 1;
                    end
                end
            % k-1 alternatyw kosumowanych
            Dsim(FirstAlts,i,j) = (MUa(FirstAlts,i,j).*(lam.^(-inva(FirstAlts)))-1).*gammas(FirstAlts);
%                 Dsim(FirstAlts,i,j) = (Lam_k_tmp.*MUa(FirstAlts,i,j) - 1).*gammas(FirstAlts); 
            end
        end
    end
end
D = mean(Dsim,3);
Err = Y-D;
end

function E = ehat(lam, MUa, inva, pg)
    E = sum(pg.*(MUa.*(lam.^(-inva))-1),1); 
end

function lam = bisect(Lam_low, Lam_up, I, error, MUa, inva, pg)
    Lam_hat = (Lam_up + Lam_low)/2;
    E_hat = ehat(Lam_hat, MUa, inva, pg);

    if abs(Lam_up - Lam_low) < error || abs(E_hat - I) < error
        lam = Lam_hat;
    else
        condb = 0;
        while condb == 0
            if E_hat < I
                Lam_up = (Lam_up + Lam_low)/2;
            else
                Lam_low = (Lam_up + Lam_low)/2;
            end
            Lam_hat = (Lam_up + Lam_low)/2;
            E_hat = ehat(Lam_hat, MUa, inva, pg);
            if abs(Lam_up - Lam_low) < error || abs(E_hat - I) < error
                condb = 1;
                lam = Lam_hat;
            end
        end
    end
end