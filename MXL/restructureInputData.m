function INPUT = restructureInputData(INPUT,EstimOpt)

INPUT.XXa = reshape(INPUT.Xa,[EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP,EstimOpt.NVarA]);
INPUT.XXa = permute(INPUT.XXa,[1 3 2]);
INPUT.YY = reshape(INPUT.Y,[EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP]);

% idx = sum(reshape(INPUT.MissingInd,[EstimOpt.NAlt,EstimOpt.NCT,EstimOpt.NP])) == EstimOpt.NAlt;
% INPUT.YYY(idx(ones(EstimOpt.NAlt,1),:,:)) = NaN; % replace YYY in missing choice-tasks with NaN
% INPUT.YY = reshape(INPUT.YYY,EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP)==1;
%INPUT.YY = reshape(INPUT.YYY,EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP);

INPUT.XXm = reshape(INPUT.Xm',[EstimOpt.NVarM,EstimOpt.NAlt*EstimOpt.NCT,EstimOpt.NP]);
% INPUT.XXm = squeeze(INPUT.XXm(:,1,:));
% if EstimOpt.NVarM == 1
%     INPUT.XXm = INPUT.XXm';
% end
INPUT.XXm = reshape(INPUT.XXm(:,1,:),[EstimOpt.NVarM,EstimOpt.NP]);

end

