function err_mtx = generateRandomDraws(EstimOpt)
% Function generates pseudo-random draws from standard normal distribution

if isfield(EstimOpt,'Seed1') == 1
   rng(EstimOpt.Seed1);
end
cprintf('Simulation with ');
cprintf('*blue',[num2str(EstimOpt.NRep) ' ']);

if EstimOpt.Draws == 1
    cprintf('*blue','Pseudo-random '); cprintf('draws \n');
    err_mtx = randn(EstimOpt.NP*EstimOpt.NRep,EstimOpt.NVarA);
elseif EstimOpt.Draws == 2 % LHS
    cprintf('*blue','Latin Hypercube Sampling '); cprintf('draws \n');
    err_mtx = lhsnorm(zeros((EstimOpt.NVarA)*EstimOpt.NP,1), ...
        diag(ones(EstimOpt.NVarA*EstimOpt.NP,1)),EstimOpt.NRep);
    err_mtx = reshape(err_mtx,[EstimOpt.NRep*EstimOpt.NP,EstimOpt.NVarA]);
elseif EstimOpt.Draws >= 3 % Quasi random draws
    if EstimOpt.Draws == 3
        cprintf('*blue','Halton ');
        cprintf(['draws (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n']);
        hm1 = haltonset(EstimOpt.NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
    elseif EstimOpt.Draws == 4 % apply reverse-radix scrambling
        cprintf('*blue','Halton '); cprintf(['draws with reverse radix scrambling (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n']);
        hm1 = haltonset(EstimOpt.NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
        hm1 = scramble(hm1,'RR2');
    elseif EstimOpt.Draws == 5
        cprintf('*blue','Sobol '); cprintf(['draws (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n']);
        hm1 = sobolset(EstimOpt.NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
    elseif EstimOpt.Draws == 6
        cprintf('*blue','Sobol '); cprintf(['draws with random linear scramble and random digital shift (skip = ',num2str(EstimOpt.HaltonSkip),'; leap = ',num2str(EstimOpt.HaltonLeap),') \n']);
        hm1 = sobolset(EstimOpt.NVarA,'Skip',EstimOpt.HaltonSkip,'Leap',EstimOpt.HaltonLeap);
        hm1 = scramble(hm1,'MatousekAffineOwen');
    end
    
    % hm1 is a representation of Unif(0, 1) -- err_mtx consists of draws from Unif(0, 1)
    err_mtx = net(hm1,EstimOpt.NP*EstimOpt.NRep); % this takes every point:
    clear hm1;
    %     err_mtx = err_mtx(:,2:EstimOpt.NVarA+1);
    
    % Transform err_mtx to draws from standard normal
    if EstimOpt.NP*EstimOpt.NRep < 3e+7
        err_mtx = icdf('Normal',err_mtx,0,1); %to be cut down later
    else % this is for very large number of draws * variables
        for i = 1:EstimOpt.NVarA
            err_mtx(:,i) = icdf('Normal',err_mtx(:,i),0,1); %to be cut down later
        end
    end
end

end

