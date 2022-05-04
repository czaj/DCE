function distribution = distType(input, fixed, input_length, type)

if nargin == 4 && strcmp(type,'lml')
    distribution = cell(length(input),1);
    for i = 1:length(input)
        switch input(i)
            case 0
                distribution(i,1) = {'a n'};
            case -1
                distribution(i,1) = {'f'};
            case 1
                distribution(i,1) = {'a l'};
            case 2
                distribution(i,1) = {'LP (n)'};
            case 3
                distribution(i,1) = {'LP (ln)'};
            case 4
                distribution(i,1) = {'Sf'};
            case 5
                distribution(i,1) = {'LSp'};
            case 6
                distribution(i,1) = {'CSp'};
            case 7
                distribution(i,1) = {'Pw CSp'};
            case 8
                distribution(i,1) = {'Pw CHISp'};
            %%%% WZ19. %%%%
            case 9
                distribution(i,1) = {'Sf'}; % Uni-Log (Reciprocal)
            case 10
                distribution(i,1) = {'Sf'}; % Pareto
            case 11
                distribution(i,1) = {'Sf'}; % Lomax
            case 12
                distribution(i,1) = {'Sf'}; % Logistic
            case 13
                distribution(i,1) = {'Sf'}; % Log-Logistic
            case 14
                distribution(i,1) = {'Sf'}; % Gumbel                  
            case 15
                distribution(i,1) = {'Sf'}; % Cauchy
            case 16
                distribution(i,1) = {'Sf'}; % Rayleigh
            case 17
                distribution(i,1) = {'Sf'}; % Exponential
            %%%% WZ19. koniec %%%%                
        end
    end
elseif nargin >= 3 && fixed == 1
    distribution = cell(input_length,1);
else
    distribution = cell(length(input),1);
    for i = 1:length(input)
        switch input(i)
            case 0
                distribution(i,1) = {'n'};
            case -1
                distribution(i,1) = {'f'};
            case 1
                distribution(i,1) = {'l'};
            case 2
                distribution(i,1) = {'S'};
            case 3
                distribution(i,1) = {'t'};
            case 4
                distribution(i,1) = {'W'};
            case 5
                distribution(i,1) = {'s-a'};
            case 6
                distribution(i,1) = {'JSb'};
            case 7
                distribution(i,1) = {'JSu'};
            %%%% WZ20. %%%%
            case 8
                distribution(i,1) = {'nic'}; % puste - nwm dlaczego (wczesniej jest 7 przypadkow)
            case 9
                distribution(i,1) = {'UL'}; % Uni-Log
            case 10
                distribution(i,1) = {'Par'}; % Pareto
            case 11
                distribution(i,1) = {'Lom'}; % Lomax
            case 12
                distribution(i,1) = {'Log'}; % Logistic
            case 13
                distribution(i,1) = {'LL'}; % Log-Logistic
            case 14
                distribution(i,1) = {'Gum'}; % Gumbel                 
            case 15
                distribution(i,1) = {'Cau'}; % Cauchy
            case 16
                distribution(i,1) = {'Ray'}; % Rayleigh
            case 17
                distribution(i,1) = {'Exp'}; % Exponential
            %%%% WZ20. koniec %%%%
        end
    end
end

end
