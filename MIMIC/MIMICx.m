function Results = MIMICx(Y, Xa, MissingInd, Xmea, Xstr, err, EstimOpt, OptimOpt)

disp(num2str(EstimOpt.Latent,'Estimating two stage HMNL with %1.0f Latent Variable(s) for starting values'))

%% Estimation of MIMIC model

tic;
EstimOpt.BActivex = EstimOpt.BActive;
EstimOpt.BActive = EstimOpt.BActivex(EstimOpt.NVarA*(1+EstimOpt.Latent)+1:end);
b0 = 0.001*ones(EstimOpt.NVarStr*EstimOpt.Latent + EstimOpt.NVarMea + EstimOpt.NVarCut,1);
LLfun = @(B) LL_mimic_MATlike(Xstr,Xmea,err,EstimOpt,OptimOpt,B);
disp(num2str(EstimOpt.Latent,'Estimating MIMIC stage with %1.0f Latent Variable(s) for starting values'))
if EstimOpt.HessEstFix == 0
	[Results.b,LL,Results.exitf,Results.output,Results.g,Results.hess] = fminunc(LLfun,b0,OptimOpt);
else
	[Results.b,LL,Results.exitf,Results.output,Results.g] = fminunc(LLfun,b0,OptimOpt);
end     

Results.LL = -LL;

%% MIMIC output 

if EstimOpt.HessEstFix == 1
	f = LL_mimic(Xstr,Xmea,err,EstimOpt,Results.b);
	Results.jacobian = numdiff(@(B) LL_mimic(Xstr,Xmea,err,EstimOpt,B),f,Results.b,isequal(OptimOpt.FinDiffType,'central'));
elseif EstimOpt.HessEstFix == 2
	Results.jacobian = jacobianest(@(B) LL_mimic(Xstr,Xmea,err,EstimOpt,B),Results.b);
elseif EstimOpt.HessEstFix == 3
	Results.jacobian = hessian(LLfun,Results.bhat);
end
if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
    Results.hess = Results.jacobian'*Results.jacobian;
    Results.ihess = inv(Results.hess);  
    Results.std = sqrt(diag(Results.ihess));     
else    
    Results.ihess = inv(Results.hess);
    Results.std = sqrt(diag(Results.ihess)); 
end   

bstr = reshape(Results.b(1:EstimOpt.NVarStr*EstimOpt.Latent),[EstimOpt.NVarStr,EstimOpt.Latent]);
LV = Xstr*bstr;
LV = LV(:,:,ones(EstimOpt.NCT*EstimOpt.NAlt,1)); % NP x Latent x NCT*NAlt
LV = permute(LV,[3 1 2]);
LV = reshape(LV,[EstimOpt.NCT*EstimOpt.NAlt*EstimOpt.NP,EstimOpt.Latent]); 

%% Estimation of MNL

disp('Estimating MNL stage for starting values')
INPUT.Xa = Xa;

for i = 1:EstimOpt.Latent
    LVi = LV(:,i);
    INPUT.Xa = [INPUT.Xa,Xa.*LVi];
end

INPUT.Y = Y;
INPUT.MissingInd = MissingInd;
EstimOpt.DISPLAY = 0;
EstimOpt.BActive = EstimOpt.BActivex(1:EstimOpt.NVarA*(1+EstimOpt.Latent));

Results.MNL_LV = MNL(INPUT,[],EstimOpt,OptimOpt);

Results.bhat = [Results.MNL_LV.bhat;Results.b];
Results.std = [Results.MNL_LV.std;Results.std];
Results.LL = Results.LL + Results.MNL_LV.LL; 

%% Output

Results.DetailsA = [Results.bhat(1:EstimOpt.NVarA),Results.std(1:EstimOpt.NVarA),pv(Results.bhat(1:EstimOpt.NVarA),Results.std(1:EstimOpt.NVarA))];
Results.DetailsL = [Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(1+EstimOpt.Latent)),Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA*(1+EstimOpt.Latent)),pv(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA*(1+EstimOpt.Latent)),Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA*(1+EstimOpt.Latent)))];
Results.DetailsS = [Results.bhat(EstimOpt.NVarA*(1+EstimOpt.Latent)+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.Latent+EstimOpt.NVarA),Results.std(EstimOpt.NVarA*(1+EstimOpt.Latent)+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.Latent+EstimOpt.NVarA),pv(Results.bhat(EstimOpt.NVarA*(1+EstimOpt.Latent)+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.Latent+EstimOpt.NVarA),Results.std(EstimOpt.NVarA*(1+EstimOpt.Latent)+1:(EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.Latent+EstimOpt.NVarA))];
Results.DetailsM = [Results.bhat((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.Latent+EstimOpt.NVarA+1:end),Results.std((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.Latent+EstimOpt.NVarA+1:end),pv(Results.bhat((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.Latent+EstimOpt.NVarA+1:end),Results.std((EstimOpt.NVarA+EstimOpt.NVarStr)*EstimOpt.Latent+EstimOpt.NVarA+1:end))];

disp(' ')
disp('MNL attributes');
disp('var.  coef.     st.err.  p-value')
disp([num2str([(1:EstimOpt.NVarA)',Results.DetailsA(:,1)],'%1.0f %8.4f'),star_sig(Results.DetailsA(:,3)),num2str(Results.DetailsA(:,2:3),'%8.4f %8.4f')])
l = 0;
for i = 1:EstimOpt.Latent
    disp(' ')
    disp(num2str(i,'Interactions with Latent Variable %1.0f'));
    disp('var.  coef.     st.err.  p-value')
    disp([num2str([(1:EstimOpt.NVarA)',Results.DetailsL(l+1:l+EstimOpt.NVarA,1)],'%1.0f %8.4f'),star_sig(Results.DetailsL(l+1:l+EstimOpt.NVarA,3)),num2str(Results.DetailsL(l+1:l+EstimOpt.NVarA,2:3),'%8.4f %8.4f')])
    l = l + EstimOpt.NVarA;
end
l = 0;
for i = 1:EstimOpt.Latent
    disp(' ')
    disp(num2str(i,'Structural equation of Latent Variable %1.0f'));
    disp('var. coef.     st.err.  p-value')
    disp([num2str([(1:EstimOpt.NVarStr)',Results.DetailsS(l+1:l+EstimOpt.NVarStr,1)],'%1.0f %8.4f'),star_sig(Results.DetailsS(l+1:l+EstimOpt.NVarStr,3)),num2str(Results.DetailsS(l+1:l+EstimOpt.NVarStr,2:3),'%8.4f %8.4f')])
    l = l + EstimOpt.NVarStr;
end

l = 0;
for i = 1:size(Xmea,2)
    disp(' ')
    disp(num2str(i,'Measurment equation %1.0f'));
    if EstimOpt.MeaSpecMatrix(i) == 0
        disp('Estimated using OLS')
        disp('var.       coef.     st.err.  p-value')
        disp([num2str(Results.DetailsM(l+1,1),'constant %7.4f'),star_sig(Results.DetailsM(l+1,3)),num2str(Results.DetailsM(l+1,2:3),'%7.4f %8.4f')])
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            k = find(EstimOpt.MeaMatrix(:,i) == 1);
            disp([num2str([k(n),Results.DetailsM(l+1+n,1)],'LV %1.f %11.4f'),star_sig(Results.DetailsM(l+1+n,3)),num2str(Results.DetailsM(l+1+n,2:3),'%8.4f %8.4f')])
           
        end
        l = l+sum(EstimOpt.MeaMatrix(:,i))+1;
        Results.DetailsM(l+1,1:3) = [exp(Results.DetailsM(l+1,1)),Results.DetailsM(l+1,2)*exp(Results.DetailsM(l+1,1)),pv(exp(Results.DetailsM(l+1,1)),Results.DetailsM(l+1,2)*exp(Results.DetailsM(l+1,1)))];
        disp([num2str(Results.DetailsM(l+1,1),'st.dev. %8.4f'),star_sig(Results.DetailsM(l+1,3)),num2str(Results.DetailsM(l+1,2:3),'%8.4f %8.4f')])
        l = l + 1;
    elseif EstimOpt.MeaSpecMatrix(i) == 1
        disp('Estimated using MNL')
        disp('var.  coef.     st.err.  p-value')
        disp([num2str([(1:(sum(EstimOpt.MeaMatrix(:,i))*(length(unique(Xmea(:,i)))-1)))',Results.DetailsM(l+1:l+(sum(EstimOpt.MeaMatrix(:,i))*(length(unique(Xmea(:,i)))-1)),1)],'%1.0f %8.4f'),star_sig(Results.DetailsM(l+1:l+(sum(EstimOpt.MeaMatrix(:,i))*(length(unique(Xmea(:,i)))-1)),3)),num2str(Results.DetailsM(l+1:l+(sum(EstimOpt.MeaMatrix(:,i))*(length(unique(Xmea(:,i)))-1)),2:3),'%8.4f %8.4f')])
        l = l + (sum(EstimOpt.MeaMatrix(:,i))*(length(unique(Xmea(:,i))) - 1));
    elseif EstimOpt.MeaSpecMatrix(i) == 2
        disp('Estimated using Ordered Probit')
        disp('var.       coef.     st.err.  p-value')
        for n = 1:sum(EstimOpt.MeaMatrix(:,i),1)
            k = find(EstimOpt.MeaMatrix(:,i) == 1);
            disp([num2str([k(n),Results.DetailsM(l+n,1)],'LV %1.0f %11.4f'),star_sig(Results.DetailsM(l+n,3)),num2str(Results.DetailsM(l+n,2:3),'%8.4f %8.4f')])
        end  
        l = l + sum(EstimOpt.MeaMatrix(:,i));
        disp([num2str([1,Results.DetailsM(l+1,1)],'Cutoff %1.0f %7.4f'),star_sig(Results.DetailsM(l+1,3)),num2str(Results.DetailsM(l+1,2:3),'%8.4f %8.4f')])
        if length(unique(Xmea(:,i))) > 2 % if attitude is not binary 
            g = [Results.DetailsM(l+1,1) ; exp(Results.DetailsM(l+2:l+length(unique(Xmea(:,i)))-1,1))];
            for n = 2:length(unique(Xmea(:,i)))-1
                stdx = sqrt(g(1:n)'*Results.ihess(EstimOpt.Latent*EstimOpt.NVarStr+l+1:EstimOpt.Latent*EstimOpt.NVarStr+l+n,EstimOpt.Latent*EstimOpt.NVarStr+ l+1:EstimOpt.Latent*EstimOpt.NVarStr+l+n)*g(1:n));
                Results.DetailsM(l+n,1:3) = [sum(g(1:n),1),stdx,pv(sum(g(1:n),1),stdx)];
                disp([num2str([n,Results.DetailsM(l+n,1)],'Cutoff %1.0f %7.4f'),star_sig(Results.DetailsM(l+n,3)),num2str(Results.DetailsM(l+n,2:3),'%8.4f %8.4f')])
            end
        end       
        l = l + length(unique(Xmea(:,i))) - 1;
    end
end

disp(' ')
disp(['LL at convergence: ',num2str(Results.LL,'%8.4f')])
disp(' ');
Results.clocknote = clock;
Results.tocnote = toc;
[~,DayName] = weekday(now,'long');
disp(['Estimation completed on ' DayName ', ' num2str(Results.clocknote(1)) '-' sprintf('%02.0f',Results.clocknote(2)) '-' sprintf('%02.0f',Results.clocknote(3)) ' at ' sprintf('%02.0f',Results.clocknote(4)) ':' sprintf('%02.0f',Results.clocknote(5)) ':' sprintf('%02.0f',Results.clocknote(6))])
disp(['Estimation took ' num2str(Results.tocnote) ' seconds ('  num2str(floor(Results.tocnote/(60*60))) ' hours ' num2str(floor(rem(Results.tocnote,60*60)/60)) ' minutes ' num2str(rem(Results.tocnote,60)) ' seconds).']);
disp(' ');