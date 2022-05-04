function stop = outputf(x,optimvalues,state)
global B_backup 
persistent LL_backup IterTime

% save tmp1

stop = false;
if isequal(state,'init')
    disp('')
	fprintf('%6s %6s %8s %16s %17s %18s %17s %12s \n','Iter.','Eval.','dB','Step','f(x)','df(x)','Opt. Cond.','Iter. time');
elseif isequal(state,'iter')
    IterTocNote = toc(IterTime);
    if optimvalues.iteration == 0
        if ~isfield('optimvalues','stepsize') || isempty(optimvalues.stepsize)
            optimvalues.stepsize = NaN;
            optimvalues.firstorderopt = NaN; % this is to make it work with fminsearch
        end
        fprintf('%4d %6d %15.10f %15.10f %19.10f %15.10f %15.10f %9.4f\n',optimvalues.iteration,optimvalues.funccount,NaN,optimvalues.stepsize,optimvalues.fval,NaN,optimvalues.firstorderopt,IterTocNote);
        B_backup = x;
        LL_backup = optimvalues.fval;
    else
        if ~isfield(optimvalues,'stepsize') || isempty(optimvalues.stepsize) % this is to make it work with fminsearch
            optimvalues.stepsize = NaN;
            optimvalues.firstorderopt = NaN; % this is to make it work with fminsearch
        end
        dB = max(abs(x - B_backup));
        if isfield(optimvalues,'procedure') && ~isempty(optimvalues.procedure)
            fprintf('%4d %6d %15.10f %15.10f %19.10f %15.10f %15.10f %9.4f    %s\n',optimvalues.iteration,optimvalues.funccount,dB,optimvalues.stepsize,optimvalues.fval,LL_backup - optimvalues.fval,optimvalues.firstorderopt,IterTocNote,optimvalues.procedure);
        else
            fprintf('%4d %6d %15.10f %15.10f %19.10f %15.10f %15.10f %9.4f\n',optimvalues.iteration,optimvalues.funccount,dB,optimvalues.stepsize,optimvalues.fval,LL_backup - optimvalues.fval,optimvalues.firstorderopt,IterTocNote);
        end
        B_backup = x;
        LL_backup = optimvalues.fval;
    end
end

IterTime = tic;


