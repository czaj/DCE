function [Results_old, EstimOpt, b0] = setStartingValues(INPUT, Results_old, EstimOpt)

global B_backup

if EstimOpt.FullCov == 0
    if exist('B_backup','var') && ~isempty(B_backup) && ...
            size(B_backup,1) == EstimOpt.NVarA*(2+EstimOpt.NVarM) + EstimOpt.NVarS + EstimOpt.NVarNLT + 2*EstimOpt.Johnson
        b0 = B_backup(:);
        disp('Using the starting values from Backup')
    elseif isfield(Results_old,'MXL_d') && isfield(Results_old.MXL_d,'b0') % starting values provided
        Results_old.MXL_d.b0_old = Results_old.MXL_d.b0(:);
        Results_old.MXL_d = rmfield(Results_old.MXL_d,'b0');
        if length(Results_old.MXL_d.b0_old) ~= ...
                EstimOpt.NVarA*2 + EstimOpt.NVarM*EstimOpt.NVarA + EstimOpt.NVarS + EstimOpt.NVarNLT + 2*EstimOpt.Johnson
            cprintf(rgb('DarkOrange'),'WARNING: Incorrect no. of starting values or model specification \n')
            Results_old.MXL_d = rmfield(Results_old.MXL_d,'b0_old');
        else
            b0 = Results_old.MXL_d.b0_old(:);
        end
    end
    if ~exist('b0','var')
        disp('Using MNL results as starting values')
        if ~(isfield(Results_old,'MNL') && isfield(Results_old.MNL,'bhat') && ...
                length(Results_old.MNL.bhat) == (EstimOpt.NVarA*(1+EstimOpt.NVarM) + EstimOpt.NVarS + EstimOpt.NVarNLT))
            EstimOpt_tmp = EstimOpt;
            EstimOpt_tmp.Display = 0;
            OptimOpt_tmp = optimoptions('fminunc');
            OptimOpt_tmp.Algorithm = 'quasi-newton';
            OptimOpt_tmp.GradObj = 'off';
            OptimOpt_tmp.Hessian = 'off';
            OptimOpt_tmp.Display = 'off';
            OptimOpt_tmp.FunValCheck= 'off';
            OptimOpt_tmp.Diagnostics = 'off';
            Results_old.MNL = MNL(INPUT,Results_old,EstimOpt_tmp,OptimOpt_tmp);
        end
        Results_old.MNL.bhat = Results_old.MNL.bhat(:);
        % b0 = [length(EstimOpt.NVarA) first parameters; length(EstimOpt.NVarA) second
        % parameters; ?? (empty)]
        b0 = [Results_old.MNL.bhat(1:EstimOpt.NVarA);...
            max(1,abs(Results_old.MNL.bhat(1:EstimOpt.NVarA)));...
            Results_old.MNL.bhat(EstimOpt.NVarA+1:end)];
        if sum(EstimOpt.Dist == 1) > 0
            if any(b0(EstimOpt.Dist == 1) < 0)
                cprintf(rgb('DarkOrange'),'WARNING: MNL estimates of log-normally distributed parameters negative - using arbitrary starting values (this may not solve the problem - sign of the attribute may need to be reversed \n')
                %b0(EstimOpt.Dist == 1 & b0(EstimOpt.Dist == 1) < 0) = 1.01;
                b0(b0(1:size(EstimOpt.Dist,2)) < 0 & EstimOpt.Dist' == 1) = 1.01;
            end
            b0(EstimOpt.Dist == 1) = log(b0(EstimOpt.Dist == 1));
        end
        if sum(EstimOpt.Dist == -1) > 0 % Fixed
            indx = find(EstimOpt.Dist == -1);
            b0(indx+EstimOpt.NVarA) = 0;   % second parameters equal to 0 (variances as default)
        end
        if sum(EstimOpt.Dist == 3) > 0 % Triangular
            indx = find(EstimOpt.Dist == 3);
            b0([indx;indx+EstimOpt.NVarA]) = ...
                [log(b0(indx) - EstimOpt.Triang');log(b0(indx) - EstimOpt.Triang')];
        end
        if sum(EstimOpt.Dist == 4) > 0 % Weibull
            indx = find(EstimOpt.Dist == 4);
            % first parameters = log(first parameters)
            % second parameters = 0
            b0([indx;indx+EstimOpt.NVarA]) = [log(b0(indx));zeros(length(indx),1)];
        end
        if sum(EstimOpt.Dist >= 5) > 0 % Johnson
            indx = find(EstimOpt.Dist >= 5);
            tmp = [b0(indx);log(b0(indx+EstimOpt.NVarA))];
            b0([indx;indx+EstimOpt.NVarA]) = [zeros(length(indx),1),ones(length(indx),1)];
            b0 = [b0;tmp];
        end
        %         else
        %             error('No starting values available - run MNL first')
    end
    
else % EstimOpt.FullCov == 1
    
    if exist('B_backup','var') && ~isempty(B_backup) && ...
            size(B_backup,1) == EstimOpt.NVarA*(1+EstimOpt.NVarM) + sum(1:EstimOpt.NVarA) + EstimOpt.NVarS + EstimOpt.NVarNLT + 2*EstimOpt.Johnson
        b0 = B_backup(:);
        disp('Using the starting values from Backup')
    elseif isfield(Results_old,'MXL') && isfield(Results_old.MXL,'b0') % starting values provided
        Results_old.MXL.b0_old = Results_old.MXL.b0(:);
        Results_old.MXL = rmfield(Results_old.MXL,'b0');
        if length(Results_old.MXL.b0_old) ~= ...
                EstimOpt.NVarA*(1+EstimOpt.NVarM) + sum(1:EstimOpt.NVarA) + EstimOpt.NVarS + EstimOpt.NVarNLT + 2*EstimOpt.Johnson
            cprintf(rgb('DarkOrange'),'WARNING: Incorrect no. of starting values or model specification \n')
            Results_old.MXL = rmfield(Results_old.MXL,'b0_old');
        else
            b0 = Results_old.MXL.b0_old;
        end
    end
    if ~exist('b0','var')
        if isfield(Results_old,'MXL_d') && isfield(Results_old.MXL_d,'bhat') && ...
                length(Results_old.MXL_d.bhat) == ((2+EstimOpt.NVarM)*EstimOpt.NVarA + EstimOpt.NVarS + EstimOpt.NVarNLT + 2*EstimOpt.Johnson)
            disp('Using MXL_d results as starting values')
            Results_old.MXL_d.bhat = Results_old.MXL_d.bhat(:);
            if sum(EstimOpt.Dist >= 3) > 0
                vc_tmp = Results_old.MXL_d.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2);
                vc_tmp(EstimOpt.Dist < 3) = vc_tmp(EstimOpt.Dist < 3).^2;
                vc_tmp = diag(vc_tmp);
            else
                vc_tmp = (diag(Results_old.MXL_d.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*2))).^2;
            end
            b0 = [Results_old.MXL_d.bhat(1:EstimOpt.NVarA);...
                vc_tmp(tril(ones(size(vc_tmp))) == 1);...
                Results_old.MXL_d.bhat(EstimOpt.NVarA*2+1:end)];
        else
            disp('Using MNL results as starting values')
            if ~(isfield(Results_old,'MNL') && isfield(Results_old.MNL,'bhat') && ...
                    length(Results_old.MNL.bhat) == (EstimOpt.NVarA*(1+EstimOpt.NVarM) + EstimOpt.NVarS + EstimOpt.NVarNLT))
                EstimOpt_tmp = EstimOpt;
                EstimOpt_tmp.Display = 0;
                OptimOpt_tmp = optimoptions('fminunc');
                OptimOpt_tmp.Algorithm = 'quasi-newton';
                OptimOpt_tmp.GradObj = 'off';
                OptimOpt_tmp.Hessian = 'off';
                OptimOpt_tmp.Display = 'off';
                OptimOpt_tmp.FunValCheck= 'off';
                OptimOpt_tmp.Diagnostics = 'off';
                Results_old.MNL = MNL(INPUT,Results_old,EstimOpt_tmp,OptimOpt_tmp);
            end
            Results_old.MNL.bhat = Results_old.MNL.bhat(:);
            b0 = [Results_old.MNL.bhat(1:EstimOpt.NVarA);...
                zeros(sum(1:EstimOpt.NVarA),1);...
                Results_old.MNL.bhat(EstimOpt.NVarA+1:end)];
            if sum(EstimOpt.Dist == 1) > 0
                if any(b0(EstimOpt.Dist == 1) < 0)
                    cprintf(rgb('DarkOrange'),'WARNING: MNL estimates of log-normally distributed parameters negative - using arbitrary starting values (this may not solve the problem - sign of the attribute may need to be reversed \n')
                    %b0(EstimOpt.Dist == 1 & b0(EstimOpt.Dist == 1) < 0) = 1.01;
                    b0(b0(1:size(EstimOpt.Dist,2)) < 0 & EstimOpt.Dist' == 1) = 1.01;
                end
                b0(EstimOpt.Dist == 1) = log(b0(EstimOpt.Dist == 1));
            end
            if sum(EstimOpt.Dist == -1) > 0 % Fixed
                b0(EstimOpt.Dist == -1) = 0;
            end
            if sum(EstimOpt.Dist == 3) > 0 % Triangular
                b0(EstimOpt.Dist == 3) = log(b0(EstimOpt.Dist == 3) - EstimOpt.Triang');
            end
            if sum(EstimOpt.Dist == 4) > 0 % Weibull
                b0(EstimOpt.Dist == 4) = log(b0(EstimOpt.Dist == 4));
            end
            if sum(EstimOpt.Dist >= 5) > 0 % Johnson
                indx = find(EstimOpt.Dist >= 5);
                tmp = b0(indx);
                b0(indx) = zeros(length(indx),1);
                b0 = [b0;tmp;zeros(length(indx),1)];
            end
            %         else
            %             error('No starting values available')
        end
    end
end


end

