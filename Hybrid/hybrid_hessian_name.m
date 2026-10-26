function name = hybrid_hessian_name(options)
% Both fminunc and fmincon expose this property; legacy HessUpdate is fminunc-only.
name = options.HessianApproximation;
if iscell(name), name = sprintf('%s(%d)',name{1},name{2}); else, name = char(name); end
end
