function stop = HMNL_outputfcn(x,optimValues,state,EstimOpt,oldOutputFcn)
% HMNL_outputfcn preserves the user's optimizer OutputFcn and optionally
% adds HMNL-specific stopping/debug messages.
%
% Important: this wrapper calls oldOutputFcn first, so DataCleanDCE's
% OptimOpt.OutputFcn = @outputf continues to print the usual iteration table.

stop = false;

% Preserve user-supplied OutputFcn output and stop flag.
if nargin >= 5 && ~isempty(oldOutputFcn)
    try
        if isa(oldOutputFcn,'function_handle')
            stop = oldOutputFcn(x,optimValues,state);
        elseif iscell(oldOutputFcn)
            for k = 1:numel(oldOutputFcn)
                if ~isempty(oldOutputFcn{k})
                    stop = oldOutputFcn{k}(x,optimValues,state) || stop;
                end
            end
        end
    catch ME
        warning('HMNL:OutputFcn','User OutputFcn failed: %s',ME.message);
    end
end

if stop || nargin < 4 || ~isstruct(EstimOpt)
    return
end

if ~isfield(EstimOpt,'MaxAbsB') || isempty(EstimOpt.MaxAbsB)
    EstimOpt.MaxAbsB = 1000;
end
if ~isfield(EstimOpt,'StopOnExtremeB') || isempty(EstimOpt.StopOnExtremeB)
    EstimOpt.StopOnExtremeB = 1;
end
% Default is 0 because outputf already prints the iteration table.
if ~isfield(EstimOpt,'ProgressOutputFcn') || isempty(EstimOpt.ProgressOutputFcn)
    EstimOpt.ProgressOutputFcn = 0;
end

switch state
    case 'init'
        if EstimOpt.ProgressOutputFcn ~= 0
            disp('HMNL optimization wrapper started. User OutputFcn remains active.');
        end
    case 'iter'
        maxAbsB = max(abs(x(:)));
        if EstimOpt.ProgressOutputFcn ~= 0
            msg = sprintf('HMNL wrapper: iter %d, eval %d, f=%.8g, max|B|=%.4g', ...
                getfield_safe(optimValues,'iteration',NaN), ...
                getfield_safe(optimValues,'funccount',NaN), ...
                getfield_safe(optimValues,'fval',NaN), ...
                maxAbsB);
            disp(msg);
        end
        if EstimOpt.StopOnExtremeB ~= 0 && maxAbsB > EstimOpt.MaxAbsB
            stop = true;
            try
                msg = sprintf('WARNING: Stopping HMNL optimization because max(abs(B)) = %.6g exceeds EstimOpt.MaxAbsB = %.6g.\n', maxAbsB, EstimOpt.MaxAbsB);
                fprintf(2,'%s',msg);
            catch
                warning('HMNL:ExtremeParameter','Stopping HMNL optimization because max(abs(B)) exceeded EstimOpt.MaxAbsB.');
            end
        end
    case 'done'
        if EstimOpt.ProgressOutputFcn ~= 0
            disp('HMNL optimization wrapper finished.');
        end
end
end

function v = getfield_safe(s,field,defaultValue)
if isstruct(s) && isfield(s,field) && ~isempty(s.(field))
    v = s.(field);
else
    v = defaultValue;
end
end
