function logP = mxl_log_panel(X,chosen,B,nAlt)
% Conditional log probability of a complete choice sequence for each column of B.
U = reshape(X*B,nAlt,[],size(B,2));
U = U - max(U,[],1);
shifted = reshape(U,size(X,1),[]);
logP = sum(shifted(chosen,:),1) - reshape(sum(log(sum(exp(U),1)),2),1,[]);
end
