function stop = outputf(x,optimvalues,state)
global B_backup
% global tolB
persistent LL_backup

% save tmp1
% return

stop = false;
% tolB = 1.e-6;
if isequal(state,'init')
    disp('')
	fprintf('%6s %6s %8s %16s %17s %18s %17s\n','Iter.','Eval.','dB','Step','f(x)','df(x)','Opt. Cond.');
end
    if isequal(state,'iter')
        if optimvalues.iteration == 0
            if isempty(optimvalues.stepsize)
                optimvalues.stepsize = 0;
            end
            fprintf('%4d %6d %15.10f %15.10f %19.10f %15.10f %15.10f\n',optimvalues.iteration,optimvalues.funccount,0,optimvalues.stepsize,optimvalues.fval,0,optimvalues.firstorderopt);
%             B_backup = [B_backup, x];
            B_backup = x;
            LL_backup = optimvalues.fval;
        else
%             dB = max(abs(x - B_backup(:,end)));
            dB = max(abs(x - B_backup));
            fprintf('%4d %6d %15.10f %15.10f %19.10f %15.10f %15.10f\n',optimvalues.iteration,optimvalues.funccount,dB,optimvalues.stepsize,optimvalues.fval,LL_backup - optimvalues.fval,optimvalues.firstorderopt);
%             B_backup = [B_backup, x];
            B_backup = x;
            LL_backup = optimvalues.fval;
%             if dB < tolB 
%                 stop = true;
%                 disp(['Exiting optimalization due to change in parameters is lower than selected tolerance (',num2str(tolB),')'])
%             end
        end
    end
end