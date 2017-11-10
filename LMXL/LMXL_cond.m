function Q = LMXL_cond(YY, XXa, b_mtx, EstimOpt)

NAlt = EstimOpt.NAlt;
NCT = EstimOpt.NCT;
NP = EstimOpt.NP;
NRep = EstimOpt.NRep;
NVarA = EstimOpt.NVarA;


WTP_space = EstimOpt.WTP_space;
WTP_matrix = EstimOpt.WTP_matrix;

if WTP_space > 0
    b_mtx(1:end-WTP_space,:) = b_mtx(1:end-WTP_space,:).*b_mtx(WTP_matrix,:);
end

b_mtx = reshape(b_mtx,NVarA,NRep,NP);

YYy = YY==1;
Q = zeros(NP, NRep);
parfor n = 1:NP
    U = reshape(XXa(:,:,n)*b_mtx(:,:,n),NAlt,NCT,NRep);
    U = exp(U - max(U)); % rescale utility to avoid exploding
    U_sum = reshape(sum(U,1),NCT,NRep);
    YYy_n = YYy(:,n);
    U_selected = reshape(U(YYy_n(:,ones(NRep,1))),NCT,NRep);
    Q(n,:) = prod(U_selected ./ U_sum,1);
end

