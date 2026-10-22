function [f,g] = probs_mdcev(data,EstimOpt,b)
% Negative log likelihood per task (MDCEV) or respondent (MMDCEV).
if exist('vce_likelihood','file') ~= 2
    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'VCE'));
end
[f,g] = vce_likelihood('MDCEV',data,EstimOpt,b,nargout > 1);
end
