% This code checks the input data and specifications that the user provides in FlexibleWtp.m

% Check for positive intergers

check_ok=0;

if ceil(NP) ~= NP | NP < 1;
   disp(['NP must be a positive integer, but it is set to ' num2str(NP)]);
   disp('Program terminated.');
   return
end

if ceil(NCS) ~= NCS | NCS < 1;
   disp(['NCS must be a positive integer, but it is set to ' num2str(NCS)]);
   disp('Program terminated.');
   return
end

if ceil(NROWS) ~= NROWS | NROWS < 1;
   disp(['NROWS must be a positive integer, but it is set to ' num2str(NROWS)]);
   disp('Program terminated.');
   return
end

% Checking XMAT %
if ( size(XMAT,1) ~= NROWS)
      disp(['XMAT has ' num2str(size(XMAT,1)) ' rows']);
      disp(['but it should have NROWS= '  num2str(NROWS)   ' rows.']);
      disp('Program terminated.');
      return
end

if sum(XMAT(:,1) > NP) ~= 0
     disp(['The first column of XMAT has a value greater than NP= ' num2str(NP)]);
     disp('Program terminated.');
     return
end

if sum(XMAT(:,1) < 1) ~= 0
     disp('The first column of XMAT has a value less than 1.');
     disp('Program terminated.');
     return
end

k=(XMAT(2:NROWS,1) ~= XMAT(1:NROWS-1,1)) & (XMAT(2:NROWS,1) ~= (XMAT(1:NROWS-1,1)+1));
if sum(k) ~= 0
    disp('The first column of XMAT does not ascend from 1 to NP.');
    disp('Program terminated.')
    return
end

if sum(XMAT(:,2) > NCS) ~= 0
     disp(['The second column of XMAT has a value greater than NCS= ' num2str(NCS)]);
     disp('Program terminated.');
     return
end

if sum(XMAT(:,2) < 1) ~= 0
     disp('The second column of XMAT has a value less than 1.');
     disp('Program terminated.');
     return
end

k=(XMAT(2:NROWS,2) ~= XMAT(1:NROWS-1,2)) & (XMAT(2:NROWS,2) ~= (XMAT(1:NROWS-1,2)+1));
if sum(k) ~= 0
    disp('The second column of XMAT does not ascend from 1 to NCS.');
    disp('Program terminated.')
    return
end


if sum(XMAT(:,3) ~= 0 & XMAT(:,3) ~= 1) ~= 0
     disp('The third column of XMAT has a value other than 1 or 0.');
     disp('Program terminated.');
     return
end

for s=1:NCS
    k=(XMAT(:,2) == s);
    if sum(XMAT(k,3)) > 1
       disp('The third column of XMAT indicates more than one chosen alternative');
       disp(['for choice situation ' num2str(s)]);
       disp('Program terminated.');
       return
    end
    if sum(XMAT(k,3)) < 1
       disp('The third column of XMAT indicates that no alternative was chosen');
       disp(['for choice situation ' num2str(s)]);
       disp('Program terminated.');
       return
    end 
end

if sum(sum(isnan(XMAT)),2) ~= 0
   disp('XMAT contains missing data.');
   disp('Program terminated.');
   return
end;

if sum(sum(isinf(XMAT)),2) ~= 0
   disp('XMAT contains an infinite value.');
   disp('Program terminated.');
   return
end;

if ceil(IDPRICE) ~= IDPRICE | IDPRICE < 4;
   disp('IDPRICE must be a positive integer >3,');
   disp('because the first three variables in XMAT cannot be explanatory variables.');
   disp(['but IDPRICE is set to ' num2str(IDPRICE)]);
   disp('Program terminated.');
   return
end

if IDPRICE > size(XMAT,2);
   disp('IDPRICE is set to a number greater than the columns in XMAT.');
   disp('Program terminated.');
   return
end

if sum(IDV(:,1) > size(XMAT,2)) ~= 0;
   disp('IDV identifies a variable that is outside XMAT.');
   disp('IDV is');
   IDV(:,1)
   disp('when each element of this vector must be no greater than')
   disp([num2str(size(XMAT,2)) ' which is the number of columns in XMAT.']);
   disp('Program terminated.');
   return
end;

if sum(IDV(:,1) <= 3) ~= 0;
   disp('Each element of IDV must exceed 3');
   disp('since the first three variables in XMAT cannot be explanatory variables.');
   disp('But IDV is');
   IDV(:,1)
   disp('which has an element below 3.')
   disp('Program terminated.');
   return
end;

if size(NAMES,2) ~= 1;
   disp(['NAMES must have 1 columns and yet it is set to have ' num2str(size(NAMES,2))]);
   disp('Be sure to separate names by semicolons.');
   disp('Program terminated.');
   return
end;

if size(NAMES,1)~=NV;
   disp(['NAMES must have the same length as IDV+1 but has length ' num2str(size(NAMES,1))]);
   disp(['while NAMES has length ' num2str(size(NAMES,1))]);
   disp('Program terminated.');
   return
end; 

if size(P_Range,2) ~= 2;
   disp(['P_Range must have 2 columns and yet it is set to have ' num2str(size(P_Range,2))]);
   disp('Program terminated.');
   return
end;

if size(P_Range,1) ~= 1;
   disp(['P_Range must have 1 row and yet it is set to have ' num2str(size(P_Range,1))]);
   disp('Program terminated.');
   return
end;

if P_Range(1,1) >= P_Range(1,2);
   disp('The second element of P_Range must exceed the first.');
   disp('Program terminated.');
   return
end;

if size(WTP_Range,2) ~= 2;
   disp(['WTP_Range must have 2 columns and yet it is set to have ' num2str(size(WTP_Range,2))]);
   disp('Program terminated.');
   return
end;

if size(WTP_Range,1) ~= size(IDV,1);
   disp('WTP_Range must have the same length as IDV.');
   disp('Program terminated.');
   return
end;

if sum(WTP_Range(1,1) >= WTP_Range(1,2))>0;
   disp('The second column of WTP_Range must exceed the first column in all rows.');
   disp('Program terminated.');
   return
end;


if ceil(NGridPts) ~= NGridPts | NGridPts < 1;
   disp(['NGridPts must be a positive integer, but it is set to ' num2str(NGridPts)]);
   disp('Program terminated.');
   return
end

if ceil(NDRAWS) ~= NDRAWS | NDRAWS < 1;
   disp(['NDRAWS must be a positive integer, but it is set to ' num2str(NDRAWS)]);
   disp('Program terminated.');
   return
end

if ceil(ThisSeed) ~= ThisSeed | ThisSeed < 1;
   disp(['ThisSeed must be a positive integer, but it is set to ' num2str(ThisSeed1)]);
   disp('Program terminated.');
   return
end

if ZTYPE ~= 1 & ZTYPE ~= 2 & ZTYPE ~= 3 ;
   disp(['ZTYPE must be 1,2, or 3, but it is set to ' num2str(ZTYPE)]);
   disp('Program terminated.');
   return
end

if ZTYPE==1
if ceil(PolyOrder) ~= PolyOrder | PolyOrder < 1;
   disp(['PolyOrder must be a positive integer, but it is set to ' num2str(PolyOrder)]);
   disp('Program terminated.');
   return
end
end

if ZTYPE==2
if ceil(NLevels) ~= NLevels | NLevels < 1;
   disp(['NLevels must be a positive integer, but it is set to ' num2str(NLevels)]);
   disp('Program terminated.');
   return
end
end

if ZTYPE==3
if ceil(NKnots) ~= NKnots | NKnots < 1;
   disp(['NKnots must be a positive integer, but it is set to ' num2str(NKnots)]);
   disp('Program terminated.');
   return
end
end

if CrossCorr ~= 0 & CrossCorr ~= 1;
   disp(['CrossCorr must be 0 or 1, but it is set to ' num2str(CrossCorr)]);
   disp('Program terminated.');
   return
end


if size(StartB,1) ~= NZ;
   disp(['StartZ must contain NZ elements in a column vector, but it has this many rows: ' num2str(size(StartB,1))]);
   disp('Program terminated.');
   return
end

if size(StartB,2) ~= 1;
   disp(['StartB must have 1 column and yet it is set to have ' num2str(size(StartB,2))]);
   disp('Be sure to separate values by semicolons.');
   disp('Program terminated.');
   return
end;

if ceil(NBins) ~= NBins | NBins < 1;
   disp(['NBins must be a positive integer, but it is set to ' num2str(NBins)]);
   disp('Program terminated.');
   return
end

if WantHessian ~= 0 & WantHessian ~= 1;
   disp(['WantHessian must be 0 or 1, but it is set to ' num2str(WantHessian)]);
   disp('Program terminated.');
   return
end

if WantBoot ~= 0 & WantBoot ~= 1;
   disp(['WantBoot must be 0 or 1, but it is set to ' num2str(WantBoot)]);
   disp('Program terminated.');
   return
end

if WantBoot==1
if ceil(NReps) ~= NReps| NReps < 1;
   disp(['NReps must be a positive integer, but it is set to ' num2str(NReps)]);
   disp('Program terminated.');
   return
end
end

if YesGPU ~= 0 & YesGPU ~= 1;
   disp(['YesGPU must be 0 or 1, but it is set to ' num2str(YesGPU)]);
   disp('Program terminated.');
   return
end

if ceil(MAXITERS) ~= MAXITERS| MAXITERS < 1;
   disp(['MAXITERS must be a positive integer, but it is set to ' num2str(MAXITERS)]);
   disp('Program terminated.');
   return
end

check_ok=1;


  


