function [data,nWorkers] = mxl_worker_data(input)
% ponytail: one resident dataset per worker; shard only if RAM limits it.
persistent previous poolPrevious workerData
pool = [];
if license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask())
    pool = gcp('nocreate');
end
if isempty(pool)
    if ~isempty(workerData) && isvalid(workerData), delete(workerData); end
    previous = [];
    poolPrevious = [];
    workerData = [];
    data = struct('Value',input);
    nWorkers = 0;
else
    if isempty(workerData) || ~isvalid(workerData) ||...
            ~isequal(pool,poolPrevious) || ~isequaln(input,previous)
        if ~isempty(workerData) && isvalid(workerData), delete(workerData); end
        workerData = parallel.pool.Constant(input);
        previous = input;
        poolPrevious = pool;
    end
    data = workerData;
    nWorkers = pool.NumWorkers;
end
end
