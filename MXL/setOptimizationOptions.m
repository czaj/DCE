function [EstimOpt, OptimOpt] = setOptimizationOptions(EstimOpt, OptimOpt, b0)

if isfield(EstimOpt,'BActive')
    EstimOpt.BActive = EstimOpt.BActive(:)';
else
    EstimOpt.BActive = ones(1,length(b0));
end

if sum(EstimOpt.Dist == -1) > 0
    if isfield(EstimOpt,'BActive') == 0 || isempty(EstimOpt.BActive)
        EstimOpt.BActive = ones(1,length(b0));
    end
    if EstimOpt.FullCov == 0
        EstimOpt.BActive(EstimOpt.NVarA+find(EstimOpt.Dist == -1)) = 0;
    elseif EstimOpt.FullCov == 1
        Vt = tril(ones(EstimOpt.NVarA));
        Vt(EstimOpt.Dist == -1,:) = 0;
        EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA)) = ...
            EstimOpt.BActive(EstimOpt.NVarA+1:EstimOpt.NVarA+sum(1:EstimOpt.NVarA)).*(Vt(tril(ones(size(Vt)))~=0)');
    end
end

if ~isempty(EstimOpt.ExpB) && EstimOpt.NumGrad == 0
    EstimOpt.NumGrad = 1;
    OptimOpt.GradObj = 'off';
    if EstimOpt.Display ~= 0
        cprintf(rgb('DarkOrange'),'WARNING: Setting user-supplied gradient off - ExpB not supported by analytical gradient \n')
    end
end

if EstimOpt.ConstVarActive == 1
    if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive) || ...
            sum(EstimOpt.BActive == 0) == 0
        error ('Are there any constraints on model parameters (EstimOpt.ConstVarActive)? Constraints not provided (EstimOpt.BActive).')
    elseif length(b0) ~= length(EstimOpt.BActive)
        error('Check no. of constraints')
    end
    disp(['Initial values: ' mat2str(b0',2)])
    disp(['Parameters with zeros are constrained to their initial values: ' mat2str(EstimOpt.BActive')])
else
    if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive) || ...
            sum(EstimOpt.BActive == 0) == 0
        EstimOpt.BActive = ones(1,length(b0));
        disp(['Initial values: ' mat2str(b0',2)])
    else
        if length(b0) ~= length(EstimOpt.BActive)
            error('Check no. of constraints')
        else
            disp(['Initial values: ' mat2str(b0',2)])
            disp(['Parameters with zeros are constrained to their initial values: ' mat2str(EstimOpt.BActive')])
        end
    end
end

end

