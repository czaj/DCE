function BActiveX = StrSelect(EstimOpt,Model)

BActive = EstimOpt.BActive(:)';
StrMatrix = EstimOpt.StrMatrix;
if Model == 0 % MIMIC
   if length(BActive) == EstimOpt.NVarstr*EstimOpt.NLatent + EstimOpt.NVarmea + EstimOpt.NVarcut
       l = 0;
       EstimOpt.NVarStr = EstimOpt.NVarstr;
   else
       error('Wrong length of EstimOpt.Bactive')
   end
elseif Model == 1 % HMNL
    if length(BActive) == (EstimOpt.NVarA*(1 + EstimOpt.NLatent + EstimOpt.NVarM) + EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut)
       l = EstimOpt.NVarA*(1 + EstimOpt.NLatent + EstimOpt.NVarM);
   else
       error('Wrong length of EstimOpt.Bactive')
   end
elseif Model == 2 % HMXL_d
    if length(BActive) == (EstimOpt.NVarA*(2 + EstimOpt.NLatent + EstimOpt.NVarM) + EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut)
       l = EstimOpt.NVarA*(2 + EstimOpt.NLatent + EstimOpt.NVarM);
   else
       error('Wrong length of EstimOpt.Bactive')
   end
elseif Model == 3 % HMXL
    if length(BActive) == (EstimOpt.NVarA*(1 + EstimOpt.NLatent + EstimOpt.NVarM) + sum(1:EstimOpt.NVarA) + EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut)
       l = EstimOpt.NVarA*(1 + EstimOpt.NLatent + EstimOpt.NVarM) + sum(1:EstimOpt.NVarA);
   else
       error('Wrong length of EstimOpt.Bactive')
   end
elseif Model == 4 % HMXL_d
    if length(BActive) == (EstimOpt.NVarA*(1 + EstimOpt.NLatent + EstimOpt.NVarM) + sum(1:(EstimOpt.NVarA + EstimOpt.NLatent)) - EstimOpt.NLatent + EstimOpt.NVarStr*EstimOpt.NLatent + EstimOpt.NVarMea + EstimOpt.NVarcut)
       l = EstimOpt.NVarA*(1 + EstimOpt.NLatent + EstimOpt.NVarM) + sum(1:EstimOpt.NVarA + EstimOpt.NLatent);
   else
       error('Wrong length of EstimOpt.Bactive')
   end
end
for i = 1:EstimOpt.NLatent
   BAcTmp = BActive(l+(i-1)*EstimOpt.NVarStr+1:l+i*EstimOpt.NVarStr);
   BAcTmp(StrMatrix(i,:) == 0) = 0;
   BActive(l+(i-1)*EstimOpt.NVarStr+1:l+i*EstimOpt.NVarStr) = BAcTmp;
end
BActiveX = BActive;