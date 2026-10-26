function [b,f,exitflag,output,g,H] = hybrid_optimize(fun,b0,options,active,startHessian)
% Optimize free coordinates, optionally whitened by a previous observed Hessian.
% The objective, user callbacks and returned derivatives use ORIGINAL parameters.
% startHessian is curvature, never a covariance matrix; final inference is separate.
b0 = b0(:); n = numel(b0);
assert(all(isfinite(b0)) && isreal(b0),'DCE:Hybrid:Start','Starting values must be finite and real.');
assert(numel(active) == n && all(ismember(active(:),[0 1])),...
    'DCE:Hybrid:Active','BActive must contain one zero/one per parameter.');
active = logical(active(:)); q = nnz(active);
assert(q > 0,'DCE:Hybrid:Active','At least one parameter must be free.');
assert(isa(options,'optim.options.Fminunc'),'DCE:Hybrid:Optimizer',...
    'Active optimization requires fminunc options (ConstVarActive=0).');
assert(isempty(options.HessianMultiplyFcn),'DCE:Hybrid:Optimizer',...
    'HessianMultiplyFcn is not supported by active-coordinate optimization.');
T = eye(q);
if nargin >= 5 && ~isempty(startHessian)
    if isequal(size(startHessian),[n n]), startHessian = startHessian(active,active); end
    assert(isequal(size(startHessian),[q q]) && isreal(startHessian) && all(isfinite(startHessian(:))),...
        'DCE:Hybrid:StartHessian','Start Hessian must match the free or full parameter vector.');
    assert(norm(startHessian-startHessian','fro') <= 1e-10*max(1,norm(startHessian,'fro')),...
        'DCE:Hybrid:StartHessian','Start Hessian must be symmetric.');
    [U,flag] = chol((startHessian+startHessian')/2);
    assert(flag == 0 && rcond(U) > 1e-10,'DCE:Hybrid:StartHessian',...
        'Start Hessian must be positive definite and well conditioned; no automatic repair.');
    T = U\eye(q);
end
oldOutput = options.OutputFcn; oldPlot = options.PlotFcn;
if ~isempty(oldOutput), options.OutputFcn = @(z,v,s) callback(oldOutput,z,v,s); end
if ~isempty(oldPlot), options.PlotFcn = @(z,v,s) callback(oldPlot,z,v,s); end
% Finite-difference settings describe solver coordinates, not original units.
for field = {'TypicalX','FiniteDifferenceStepSize'}
    value = options.(field{1});
    if isnumeric(value) && numel(value) > 1
        assert(numel(value) == n,'DCE:Hybrid:DifferenceOptions',...
            '%s must be scalar or full length.',field{1});
        options.(field{1}) = value(active);
    end
end
needHessian = nargout >= 6;
if needHessian
    [z,f,exitflag,output,gz,Hz] = fminunc(@objective,zeros(q,1),options);
else
    [z,f,exitflag,output,gz] = fminunc(@objective,zeros(q,1),options);
end
b = expand(z); g = zeros(n,1); g(active) = T'\gz;
H = [];
if needHessian && ~isempty(Hz)
    H = zeros(n); H(active,active) = (T'\Hz)/T;
end
output.hybridActiveParameters = q;
output.hybridPreconditioned = nargin >= 5 && ~isempty(startHessian);
output.hybridCoordinateFirstOrderOpt = output.firstorderopt;
output.firstorderopt = max(abs(g(active)));

    function x = expand(z)
        x = b0; x(active) = b0(active)+T*z;
    end
    function [value,gradient,hessian] = objective(z)
        if nargout >= 3
            [value,gg,hh] = fun(expand(z));
            hessian = T'*hh(active,active)*T;
        elseif nargout >= 2
            [value,gg] = fun(expand(z));
        else
            value = fun(expand(z));
        end
        if nargout >= 2, gradient = T'*gg(active); end
    end
    function stop = callback(functions,z,values,state)
        if ~iscell(functions), functions = {functions}; end
        if isfield(values,'gradient') && numel(values.gradient) == q
            originalGradient = zeros(n,1); originalGradient(active) = T'\values.gradient(:);
            values.gradient = originalGradient;
        end
        if isfield(values,'searchdirection') && numel(values.searchdirection) == q
            direction = zeros(n,1); direction(active) = T*values.searchdirection(:);
            values.searchdirection = direction;
        end
        stop = false;
        for k = 1:numel(functions)
            if ~isempty(functions{k}), stop = functions{k}(expand(z),values,state) || stop; end
        end
    end
end
