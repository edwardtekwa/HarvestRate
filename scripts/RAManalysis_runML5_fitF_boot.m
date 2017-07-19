%Edward Tekwa 2/19/17, 3/23/17, 5/24/17
%run maximum likelihood procedure to fit harvest rate evolution model (with
%imperfection) to RAM database fisheries' changes in harvest rates
%first, run RAManalysis_boot_pooled5_minTrans.m, or load:
%load RAManalysis_runML5_workspace.mat
load ('RAManalysis fit 10000 F permutations xi model.mat') %this contains both original data for analysis and final analysis outputs for reference


%fit model to harvest rates (1), changes in harvest rate (2), or both (3)
fitType=1;

%fit to raw data (1), F<2 (2), log2(F) (3), or log2(F)<1 (4)
dataType=3;
%limit fitting data to F<2
% avgAllPt2=avgAllPt(:,:,avgAllPt(:,1,:)<2);
% avgAllChange2=avgAllChange(:,:,avgAllPt(:,1,:)<2);

avgAllChange2=avgAllChange;
%avgAllChange2=reshape((diff(XY_Yr(:,1,:))),1,1,[]);
avgAllChange2=avgAllChange2(~isnan(avgAllChange2(1,1,:)));
eFchange=avgAllChange2(:,1,:);
firstPt=[];
avgAllPt2=[];
Fishery_Yr2={};
j=1;
for i=1:length(XY_Yr)
    %posFirstPt=min(find(~isnan(XY_Yr(:,1,i)))); %first known time point of current fishery
    posFirstPt=min(find(~isnan(XY_Yr(:,1,i)).*XY_Yr(:,1,i)>0));
    if nanmean(XY_Yr(posFirstPt:end,5,i))>0
        firstPt(j,:)=XY_Yr(posFirstPt,1,i);
        avgAllPt2(:,:,j)=nanmean(XY_Yr(posFirstPt:end,:,i));
        Fishery_Yr2{j,1}=Fishery_Yr{i,1};
        Fishery_Yr2{j,2}=Fishery_Yr{i,2};
        j=j+1;
    end
end
%avgAllPt2=avgAllPt; %use averages of each fishery
%avgAllPt2=reshape(XY(:,1:2)',1,2,[]); %use all points
validPos=avgAllPt2(1,1,:)~=0 & log(avgAllPt2(1,5,:))>0; %find fisheries where harvest rate is not zero and ln(MSY)>0
avgAllPt2=avgAllPt2(1,:,validPos);
firstPt=firstPt(validPos,1);
Fishery_Yr2=Fishery_Yr2(validPos,:);
if dataType==3
    firstPt=log2(firstPt);
    avgAllPt2(:,1,:)=log2(avgAllPt2(:,1,:)); %use log2(F)
end
L=avgAllPt2(:,2,:);
eF=avgAllPt2(:,1,:);
allMSY=avgAllPt2(:,5,:);

%plot harvest rate F as a function of cost/benefit and MSY
%obtain linear regression
x1=reshape(L,[],1)/std(L);
%x2=reshape(log(allMSY+1),[],1);
x2=reshape(log(allMSY*1000),[],1)/std(log(allMSY*1000));
x3=reshape(firstPt>0,[],1)/std(firstPt>0); %binary: 1 if F_init >F_msy, 0 otherwise
%X=[ones(size(x1)) x1 x2 x1.*x2];
%X=[ones(size(x1)) x1 x2];
X=[ones(size(x1)) x1 x2 x3];
[b, bint, r, rint, linstats]=regress(reshape(eF,[],1),X);
LinRegressCoeffs=b %intercept, I/V coefficient, MSY coefficient, init F coefficient
LinRegressIntervals=bint
linstats %R^2, F, p, err var

figure('Color', [1 1 1]);
scatter(firstPt>0,eF);
hold on
refline(1,0);
refhor=refline(0,0);
set(refhor,'Color','k');
plot([0 0],ylim,'k')
xlabel 'initial log_2(F/F_{MSY})'
ylabel 'mean log_2(F/F_{MSY}) (excludes first year)'

scrsz = get(0,'ScreenSize');
figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3) scrsz(4)/2]);
subplot(1,3,3)
histogram(avgAllPt2(:,1,:));
ylabel 'frequency'
xlabel 'harvest rate log_2(F/F_{MSY})'
subplot(1,3,2)
histogram(log(allMSY*1000));
ylabel 'frequency'
xlabel 'fish stock ln(MSY[kg]))'
subplot(1,3,1)
histogram(avgAllPt2(:,2,:));
ylabel 'frequency'
xlabel 'cost/benefit \gamma'

%freeparstart=[0.0222,0.01,3.9543,1]; %[E xi R b]=[imperfection term, evol rate, r*Smsy, beta(risk aversion)]
%freeparstart=[2.61,1,4.9307,1];
%freeparstart=[2.7897    1.0000    4.5784    1.0000]; %fval=9.9474
%freeparstart=[0.00435070286453954,1,2.75676703571226,1]; %fval=20.5125 (fit to harvest rates<2, weight exp=0)
%freeparstart=[0.0175386828353086,1,2.64315077123759,1]; %fval=0.8214
%freeparstart=[0.370643758745823,0]; %delta, imperfection
%freeparstart=[0.0757902804248014,0];
%freeparstart=[0.0757903737291numPts,1000];
freeparstart=[10,1000]; %xi, MSY unit scaling [/tonnes]
freeparmin=[0 1000];
freeparmax=[1000 1000]; %[100 1000]
% freeparstart=[1,100];
% freeparmin=[1 0];
% freeparmax=[1 1e6]; %[100 1000]
freeparmin2=[0 0];
freeparmax2=[100 1e6];

fval=100000;
eF=avgAllPt2(:,1,:); %mean harvest rate over years
eF2=avgAllPt2(:,1,:);
L=avgAllPt2(:,2,:); %mean (cost+subsidy)/landing value over years
L2=avgAllPt2(:,2,:);
%V=avgAllPt2(:,3,:); %mean landing value/landings over years
Imin=min(L);
Imax=max(L);
%construct weight of each data point as the inverse of number of points in
%its F and L window
WeightExp=0; %exponent determines relative weight to each datapoint (0=equal for all, 1=inverse to # nearby datapoints (weigh up rare points), -1=proportional to # nearby datapoints (weigh up abundant points), etc.)
PtWeight=ones(1,length(L));
maxPtWeight_to_total=max(PtWeight)/sum(PtWeight)

PtWeight2=ones(1,length(avgAllChange2));
maxPtWeight_to_total2=max(PtWeight2)/sum(PtWeight2)

%calculate null model statistics
MeanChange=nanmean(avgAllChange2(:,1,:))
%WeightedMeanChange=nansum(PtWeight2.*reshape(avgAllChange2(:,1,:),1,[]))/sum(PtWeight)
STotDirNoWeightMean=nansum(reshape((nanmean(avgAllChange2(:,1,:)).*avgAllChange2(:,1,:))<=0,1,[])) %total weighted sum of squares (with weighted directional differences of null model) 
%STotDirMean=nansum(PtWeight2.*reshape((nanmean(avgAllChange2(:,1,:)).*avgAllChange2(:,1,:))<=0,1,[])) %total weighted sum of squares (with weighted directional differences of null model) 
STotDirMean=10000;
%STotDir2Mean=nansum(PtWeight2.*reshape((WeightedMeanChange.*avgAllChange2(:,1,:))<=0,1,[])) %total weighted sum of squares (with weighted directional differences of null model) 
%STotMag=nansum(PtWeight2.*reshape(((nanmean(avgAllChange2(:,1,:))-avgAllChange2(:,1,:)).^2),1,[])) %total weighted sum of squares (with weighted magnitude differences of null model)
%STotMag2=nansum(PtWeight2.*reshape(((WeightedMeanChange-avgAllChange2(:,1,:)).^2),1,[])) %total weighted sum of squares (with weighted magnitude differences of null model)
STotMag=10000;
STotmeans=nansum(reshape((nanmean(eF)-eF).^2,1,[])) %total weighted sum of squares
%STotmeans=sum(~isnan(eF))*nanvar(eF); %total sum of squares
%LLweightedmeans=sum(PtWeight.*LLmeans);
Mean=nanmean(eF)
SD=std(eF)
WeightedMean=nansum(PtWeight.*reshape(avgAllPt2(:,1,:),1,[]))/sum(PtWeight)
SWeightedmeans=nansum(PtWeight.*reshape((WeightedMean-eF).^2,1,[])) %total weighted sum of squares
LLInit=nansum(reshape((firstPt(:,1)-eF(:)).^2,1,[])) %total weighted sum of squares of initial condition model
InitCondR2=1-LLInit/SWeightedmeans

%evaluate separate line model (instead of using means of F and dF, use
%constant models that fit F and dF separately, for current weight
%[ParamEstsdF,fvaldF,exitflag,output] = fminsearchbnd(@(params) RAManalysis_ML2(avgAllPt2,avgAllChange2,PtWeight2,Window,Fwin,params),freeparstart2,freeparmin0,freeparmax0,options);
STotDir=10000;

global prevLLdF
global prevLLF
prevLLdF=STotDir;
prevLLF=SWeightedmeans;

%evaluate null model for current weight
%[ParamEstsNull,fvalNull,exitflag,output] = fminsearchbnd(@(params) RAManalysis_ML2_fitF_dF(avgAllPt2,avgAllPt2,avgAllChange2,PtWeight2,PtWeight2,Window,Fwin,1,params),freeparstart0,freeparmin0,freeparmax0,options);
%STotDirNull=RAManalysis_ML2(avgAllPt2,avgAllChange2,PtWeight2,Window,Fwin,ParamEstsNull)
%fval=(LLdF/prevLLdF)+(LLF/prevLLF), so LLF=(fval-(LLdF/prevLLdF))*prevLLF:
%STotmeansNull=RAManalysis_ML_fitF(avgAllPt2,PtWeight2,Window,Fwin,dataType,ParamEstsNull);
fvalNull=10000;
STotDirNull=10000;
STotmeansNull=(fvalNull-(STotDirNull/prevLLdF))*prevLLF %this is faster
%F0=1-1/ParamEstsNull(1)+((ParamEstsNull(1)^2+1)^(1/2))/ParamEstsNull(1); %null model's F expectation
if fitType==1
    STotmeansNull=SWeightedmeans; %if fitting only to F
    F0=WeightedMean;
elseif fitType==2
    STotDirNull=STotDir; %if fitting only to dF
end
F0
options = optimset('MaxFunEvals', 2000,'Display','off','TolFun',1e-4,'TolX',1e-4); %for fminsearch

% disp('exploring parameter space:');
% numSamples=50;
% LLmap=ones(numSamples+1,numSamples+1)*NaN;
% for i=1:numSamples+1
%     for j=1:numSamples+1
%         freeparstartN=[freeparmin(1)+(freeparmax(1)-freeparmin(1))*(i-1)/numSamples freeparmin(2)+(freeparmax(2)-freeparmin(2))*(j-1)/numSamples]; %[E xi R b]=[imperfection term, evol rate, r*Smsy, beta(risk aversion)]
%         %fvalN=RAManalysis_ML_fitF(avgAllPt2,firstPt,PtWeight,Window,Fwin,F0,dataType,freeparstartN);
%         fvalN =RAManalysis_ML_fitF_msy(avgAllPt2,firstPt,freeparstartN);
%         LLmap(i,j)=1-fvalN/STotmeansNull;
%     end
% end
% figure;
% LLmap2=LLmap;
% LLmap2(LLmap<0)=0;
% contourf(LLmap2);
% xlabel '\eta';
% ylabel 'MSY unit [1000kg/#]';
% title 'R^2';
% %caxis([min(LLmap(:)) min(LLmap(:))+200]);
% xticks([1:numSamples/10:numSamples+1]);
% xticklabels([freeparmin(1):(freeparmax(1)-freeparmin(1))/10:freeparmax(1)]);
% yticks([1:numSamples/10:numSamples+1]);
% yticklabels([freeparmin(2):(freeparmax(2)-freeparmin(2))/10:freeparmax(2)]);

figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3) scrsz(4)/2]);
%fit eta with different MSY units:
etasR2=[];
etas=[];
for i=1:100
    freeparminN=[0 10^((i-1)/10-1)]; %start at 10000kg
    freeparmaxN=[100 10^((i-1)/10-1)];
    [ParamEstsN,fvalN,exitflag,output] = fminsearchbnd(@(params) RAManalysis_ML_fitF_msy(avgAllPt2,firstPt,params),freeparstart,freeparminN,freeparmaxN,options);
    etasR2(i)=1-fvalN/STotmeansNull;
    etas(i)=ParamEstsN(1);
end
subplot(1,2,1)
yyaxis left
plot(etasR2,'LineWidth',2)
ylabel 'R^2'
yyaxis right
plot(etas,'LineWidth',2)
ylabel '\xi'
title 'model R^2 and \xi estimate at various MSY units'
xticks([1:10:100])
xticklabels(100000./(10.^([1:1:10])))
xlabel 'MSY unit in kg'

%evaluate fit of different eta values at MSY unit=[kg]
MSYsR2=[];
for i=60:160
    freeparstartN=[i/10 1000]; %start at 10000kg
    fvalN =RAManalysis_ML_fitF_msy(avgAllPt2,firstPt,freeparstartN);
    MSYsR2(i)=1-fvalN/STotmeansNull;
end
subplot(1,2,2)
plot(MSYsR2,'k','LineWidth',2)
xticks([60:20:160])
xticklabels([60:20:160]/10)
xlabel '\xi'
ylabel 'R^2'
title 'R^2 of model with various \xi values and MSY unit=[kg]'

disp('fitting full model by estimating both \xi and \MSY unit:');
[ParamEstsN,fvalN,exitflag,output] = fminsearchbnd(@(params) RAManalysis_ML_fitF_msy(avgAllPt2,firstPt,params),freeparstartN,freeparmin2,freeparmax2,options);
R2N=1-fvalN/STotmeansNull;
disp([num2str(i) ',R^2=' num2str(R2N) ', LL=' num2str(fvalN) ', \xi=' num2str(ParamEstsN(1)) ', MSY unit [1000kg/#]=' num2str(ParamEstsN(2))]);


disp('fitting full model by estimating \xi only:');
for i=1:50
    prevLLdF=STotDir;
    prevLLF=SWeightedmeans;
    freeparstartN=freeparmin+(freeparmax-freeparmin).*rand(1,2); %[E xi R b]=[imperfection term, evol rate, r*Smsy, beta(risk aversion)]
    if i==1
        freeparstartN=freeparstart;
    end
    %to fit model to harvest rates, use all data:
    [ParamEstsN,fvalN,exitflag,output] = fminsearchbnd(@(params) RAManalysis_ML_fitF_msy(avgAllPt2,firstPt,params),freeparstartN,freeparmin,freeparmax,options);
     R2N=1-fvalN/STotmeansNull;
    disp([num2str(i) ',R^2=' num2str(R2N) ', LL=' num2str(fvalN) ', \xi=' num2str(ParamEstsN(1)) ', MSY unit [1000kg/#]=' num2str(ParamEstsN(2))]);
    if fvalN<fval
        fval=fvalN;
        freeparstart=freeparstartN;
        ParamEsts=ParamEstsN;
    end
end
freeparstart
ParamEsts %=freeparstart
fval

freeparstart=ParamEsts;


%evaluate model
%LL=RAManalysis_ML2(avgAllPt2,avgAllChange2,PtWeight2,Window,Fwin,ParamEsts);
LLweightedmeans=fval;%R2ChangeModel_direction=1-LL/STotDir
%R2ChangeModel_direction_vsNull=1-LL/STotDirNull
%R2means=1-LLweightedmeans/SWeightedmeans
R2means_vsNull=1-LLweightedmeans/STotmeansNull
%R2tot

%fit model variant with only lower equilibrium (economic optimum model)
[ParamEstsN,fvalN,exitflag,output] = fminsearchbnd(@(params) RAManalysis_ML_fitF_msy_lowerOnly(avgAllPt2,firstPt,params),freeparstartN,freeparmin,freeparmax,options);
R2N=1-fvalN/STotmeansNull;
disp(['lower solution only: R^2=' num2str(R2N) ', LL=' num2str(fvalN) ', \xi=' num2str(ParamEstsN(1)) ', MSY unit [1000kg/#]=' num2str(ParamEstsN(2))]);

%fit model variant with only upper equilibrium
[ParamEstsN,fvalN,exitflag,output] = fminsearchbnd(@(params) RAManalysis_ML_fitF_msy_upperOnly(avgAllPt2,firstPt,params),freeparstartN,freeparmin,freeparmax,options);
R2N=1-fvalN/STotmeansNull;
disp(['upper solution only: R^2=' num2str(R2N) ', LL=' num2str(fvalN) ', \xi=' num2str(ParamEstsN(1)) ', MSY unit [1000kg/#]=' num2str(ParamEstsN(2))]);


%plot data and model predictions at log10(MSY+1)=[1:6], coded with marker
%color
figure; scatter3(L,log(allMSY*1000),eF,[],reshape(log(allMSY*1000),[],1),'jitter','on')
xlabel 'total cost/total landing value (\gamma)'
ylabel 'log(MSY[kg])'
zlabel 'log2(F/F_{MSY})'


%randomize data, reapply weighting according to original rolling windows,
%and assess new model fit to the randomized data (R^2 distribution)
R2RandDist=[];
R2MagRandDist=[];
R2FRandDist=[];
R2RandDist_vsNull=[];
R2FRandDist_vsNull=[];
RDist=[];
EDist=[];
display(['Model R^2: ' num2str(R2means_vsNull)]);
display('starting randomization:');
numRand=100000;
for i=1:numRand
%     randL=datasample(L2,length(L2)); %randomize I/V data
%     randFchange=datasample(eFchange,length(eFchange)); %resample with replacement mean harvest rate changes
%     randF=datasample(eF,length(eF)); %randomized harvest rates
    randL=L2; %original I/V data
    randFchange=eFchange(randperm(size(eFchange,3))); %randomized mean harvest rate changes
    randFOrder=randperm(size(eF,3));
    randF=eF(randFOrder); %randomized harvest rates
    %randFirstPt=firstPt(randperm(size(firstPt,1))); %randomized first harvest rate of each fishery
    randFirstPt=firstPt(randFOrder); %preserve initial condition with fishery
    %randMSY=allMSY(randperm(size(allMSY,3))); %randomize stock MSY
    randMSY=allMSY; %preserve MSY with I/V
    PtWeightBoot=ones(1,length(randL));
    RandSTotDir2=10000;
    WeightedRandMean=nansum(PtWeightBoot.*reshape(randF,1,[]))/sum(PtWeightBoot);
    FBoot0=WeightedRandMean;
    RandFTot=nansum(PtWeightBoot.*reshape((WeightedRandMean-randF).^2,1,[])); %total weighted sum of squares
    prevLLdF=RandSTotDir2; %update estimated variances to weigh objectives in subsequent fits within this boot
    prevLLF=RandFTot;
    STotDirBootNull=RandSTotDir2;
    STotmeansBootNull=RandFTot;

    %fit full model to randomized data
    if fitType==1
        [ParamEstsBoot,fvalBoot,exitflag,output] = fminsearchbnd(@(params) RAManalysis_ML_fitF_msy([randF randL zeros(1,1,length(randF)) zeros(1,1,length(randF)) randMSY],randFirstPt,params),freeparstart,freeparmin,freeparmax,options);
        for j=0:0 %extra checks: try random inital guesses
            freeparstartN=freeparmin+(freeparmax-freeparmin).*rand(1,2); %random initial guess
            [ParamEstsBootN,fvalBootN,exitflag,output] = fminsearchbnd(@(params) RAManalysis_ML_fitF_msy([randF randL zeros(1,1,length(randF)) zeros(1,1,length(randF)) randMSY],randFirstPt,params),freeparstartN,freeparmin,freeparmax,options);
            if fvalBootN<fvalBoot
                fvalBoot=fvalBootN;
                ParamEstsBoot=ParamEstsBootN;
            end
        end
        %randLL=Inf;
        randFLL=fvalBoot;
    elseif fitType==2
        [ParamEstsBoot,fvalBoot,exitflag,output] = fminsearchbnd(@(params) RAManalysis_ML2([randF randL],randFchange,PtWeightBoot,Window,Fwin,params),freeparstart,freeparmin,freeparmax,options);
        randLL=fvalBoot;
        randFLL=RAManalysis_ML_fitF([randF randL],randFirstPt,PtWeightBoot,Window,Fwin,FBoot0,dataType,ParamEstsBoot);   
    else
        [ParamEstsBoot,fvalBoot,exitflag,output] = fminsearchbnd(@(params) RAManalysis_ML2_fitF_dF([randF randL],[randF randL],randFchange,PtWeightBoot,PtWeightBoot,Window,Fwin,FBoot0,params),freeparstart,freeparmin,freeparmax,options);
        randLL=RAManalysis_ML2([randF randL],randFchange,PtWeightBoot,Window,Fwin,ParamEstsBoot);
        randFLL=(fvalBoot-(randLL/prevLLdF))*prevLLF; %this is faster
end
    R2RandDist(i)=1-randLL/RandSTotDir2;
    R2RandDist_vsNull(i)=1-randLL/STotDirBootNull;
    R2FRandDist(i)=1-randFLL/RandFTot;
    R2FRandDist_vsNull(i)=1-randFLL/STotmeansBootNull;
    RDist(i)=ParamEstsBoot(1);
    EDist(i)=ParamEstsBoot(2);
    display([num2str(i) ' F R^2: ' num2str(R2FRandDist_vsNull(i)) ', LL=' num2str(randFLL) ', R=' num2str(ParamEstsBoot(1)) ', E=' num2str(ParamEstsBoot(2))]);
end
%pFModel=sum(R2FRandDist>=R2means)/length(R2RandDist)
pFModel_vsNull=sum(R2FRandDist_vsNull>=R2means_vsNull)/length(R2RandDist)
pXi=2*min(sum(ParamEsts(1)>=RDist)/length(RDist),sum(ParamEsts(1)<RDist)/length(RDist)) %2-sided test
%pE=2*min(sum(ParamEsts(2)>=EDist)/length(EDist),sum(ParamEsts(2)<EDist)/length(EDist)) %2-sided test
meanNullR2=mean(R2FRandDist_vsNull)
stdNullR2=std(R2FRandDist_vsNull)
pNullR2=sum(R2FRandDist_vsNull<=0)/length(R2RandDist)
dR2dist=R2means_vsNull-R2FRandDist_vsNull;
mean_dR2=mean(dR2dist)
std_dR2=std(dR2dist)
pFModel_vsNull_minusRand=sum(dR2dist<=0)/length(dR2dist)


scrsz = get(0,'ScreenSize');
figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3) scrsz(4)/2]);
MSYlevel=[12 14 16 18 20];

%subplot(3,2,5)
subplot(1,2,1)
his=histogram(R2FRandDist_vsNull,20);
hold on
line([R2means_vsNull R2means_vsNull], [0 max(his.Values)],'LineWidth',2,'Color','k')
title({['modelled harvest rates R^2'] ['vs. randomized data fit']})
xlabel(['R^2 (=' num2str(R2means_vsNull,3) ', p=' num2str(pFModel_vsNull,3) ')'])
ylabel(['frequency (' num2str(numRand) ' resamples)'])

%subplot(3,2,6)
subplot(1,2,2)
his=histogram(RDist,20);
hold on
line([ParamEsts(1) ParamEsts(1)], [0 max(his.Values)],'LineWidth',2,'Color','k')
title({['estimated \xi'] ['vs. randomized data fit']})
xlabel(['\xi (=' num2str(ParamEsts(1),3) ', p=' num2str(pR,3) ')'])

figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/2]);
his=histogram(dR2dist,20);
title({['modelled harvest rates \DeltaR^2'] ['vs. randomized data fit']})
xlabel(['mean \DeltaR^2 (=' num2str(mean_dR2,3) ', p=' num2str(pFModel_vsNull_minusRand,3) ')'])
ylabel(['frequency (' num2str(numRand) ' resamples)'])

theta=ParamEsts(1);
V=1;

% subplot(3,2,1:4)
% colormap hsv
% Colors=colormap;
% cmapLength=max(log(allMSY*1000))-min(log(allMSY*1000));
% scatter3(L2,log(allMSY*1000),eF,[],reshape(log(allMSY*1000),[],1),'filled','jitter','on')
% hold on
% for jMSY=1:5 %construct solutions for log10(MSY+1)=[1:6]
%     MSYs=ones(1,length(Lrange))*MSYlevel(jMSY);
%     Fsols=ones(3,length(Lrange))*(NaN);
%     for i=1:length(Lrange)
%         eL=Lrange(i);
%         %Fsol=[0; log(1+(1-1/(L*delta*jMSY))^0.5)/log(2); log(1-(1-1/(L*delta*jMSY))^0.5)/log(2)]; %equation without E (faster)
%         %Fsol=[0; log(1+(1-1/(L*delta*(MSYlevel(jMSY)-log(delta))))^0.5)/log(2); log(1-(1-1/(L*delta*(MSYlevel(jMSY)-log(delta))))^0.5)/log(2)]; %equation without E (faster)
%         Fsol=[0; log(1+(1-1/(eL*(1/theta)*MSYlevel(jMSY)))^0.5)/log(2); log(1-(1-1/(eL*(1/theta)*MSYlevel(jMSY)))^0.5)/log(2)]; %equation without E (faster)
%         Fsol(imag(Fsol)~=0)=-100;
%         Fsol=sort(Fsol,'descend');
%         Fsol(Fsol==-100)=NaN;
%         Fsols(1:length(Fsol),i)=Fsol;
%     end
%     plot1=plot3(Lrange,MSYs,Fsols','Color',Colors(ceil(64*(MSYlevel(jMSY)-min(log(allMSY*1000)))/cmapLength),:,:),'LineWidth',2);
%     plot1(2).LineStyle='--';
% end
% view(78,18)
% xlim([0.6 2.25])
% zlim([min(eF2)-.05 max(eF2)+0.05])
% zlabel 'log_2(F/F_M_S_Y)'
% xlabel 'total cost/total landing value (\gamma)'
% ylabel 'ln(MSY[kg])'
% cbar=colorbar;
% set(get(cbar,'ylabel'), 'String', 'ln(MSY[kg])');
% title('fitted model')


%plot model fit at various MSYs
scrsz = get(0,'ScreenSize');
figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3) scrsz(4)]);
colormap hsv
subplot(2,3,1)
indicesBelowMSY=find(reshape(firstPt,1,1,[])<0);
indicesAboveMSY=find(reshape(firstPt,1,1,[])>=0);
numInitBelowMSY=length(indicesBelowMSY)
numInitAboveMSY=length(indicesAboveMSY)
%scatter3(L2,log(allMSY*1000),eF,[],reshape(log(allMSY*1000),[],1),'filled','jitter','on')
%scatter3(L2,log(allMSY*1000),eF,'k','filled','jitter','on')
scatter3(L2(indicesAboveMSY),log(allMSY(indicesAboveMSY)*1000),eF(indicesAboveMSY),'ok','filled','jitter','on','jitterAmount',0.02)
hold on
scatter3(L2(indicesBelowMSY),log(allMSY(indicesBelowMSY)*1000),eF(indicesBelowMSY),'ok','jitter','on','jitterAmount',0.02)
    
hold on
for jMSY=1:5 %construct solutions for log10(MSY+1)=[1:6]
    MSYs=ones(1,length(Lrange))*MSYlevel(jMSY);
    Fsols=ones(3,length(Lrange))*(NaN);
    for i=1:length(Lrange)
        eL=Lrange(i);
        %Fsol=[0; log(1+(1-1/(L*delta*jMSY))^0.5)/log(2); log(1-(1-1/(L*delta*jMSY))^0.5)/log(2)]; %equation without E (faster)
        %Fsol=[0; log(1+(1-1/(L*delta*(MSYlevel(jMSY)-log(delta))))^0.5)/log(2); log(1-(1-1/(L*delta*(MSYlevel(jMSY)-log(delta))))^0.5)/log(2)]; %equation without E (faster)
        Fsol=[0; log(1+(1-1/(eL*(1/theta)*MSYlevel(jMSY)))^0.5)/log(2); log(1-(1-1/(eL*(1/theta)*MSYlevel(jMSY)))^0.5)/log(2)]; %equation without E (faster)
        Fsol(imag(Fsol)~=0)=-100;
        Fsol=sort(Fsol,'descend');
        Fsol(Fsol==-100)=NaN;
        Fsols(1:length(Fsol),i)=Fsol;
    end
    plot1=plot3(Lrange,MSYs,Fsols','Color','k','LineWidth',2);
    plot1(2).LineStyle='--';
end
view(90,0)
xlim([0.6 2.25])
zlim([min(eF2)-.05 max(eF2)+0.05])
zlabel 'log_2(F/F_M_S_Y)'
xlabel 'cost/benefit (\gamma)'
ylabel 'ln(MSY[kg])'
title('fitted model')

%record fishery information for each MSY bin, sorted by increasing
%cost/benefit and decreasing mean F/Fmsy
%(gamma) within each bin
FisheryBins={}; %columns 1-7: stock, countries, cost/benefit, MSY, initial harvest rate, mean harvest rate, index according to Fishery_Yr2

for jMSY=1:5 %construct solutions for log10(MSY+1)=[1:6]
    subplot(2,3,jMSY+1)
    indicesBelowMSY=find(log(allMSY*1000)>MSYlevel(jMSY)-1 & log(allMSY*1000)<=MSYlevel(jMSY)+1 & reshape(firstPt,1,1,[])<0);
    indicesAboveMSY=find(log(allMSY*1000)>MSYlevel(jMSY)-1 & log(allMSY*1000)<=MSYlevel(jMSY)+1 & reshape(firstPt,1,1,[])>=0);
    FisheryBins(:,:,jMSY)={Fishery_Yr2([indicesAboveMSY;indicesBelowMSY],:)};
    sortrows(
    scatter(L2(indicesAboveMSY),eF(indicesAboveMSY),'ok','filled','jitter','on','jitterAmount',0.02)
    hold on
    text(L2([indicesAboveMSY;indicesBelowMSY]),eF([indicesAboveMSY;indicesBelowMSY]), Fishery_Yr2([indicesAboveMSY;indicesBelowMSY],1), 'horizontal','left', 'vertical','bottom','FontSize',5)
    scatter(L2(indicesBelowMSY),eF(indicesBelowMSY),'ok','jitter','on','jitterAmount',0.02)
    MSYs=ones(1,length(Lrange))*MSYlevel(jMSY);
    Fsols=ones(3,length(Lrange))*(NaN);
    for i=1:length(Lrange)
        eL=Lrange(i);
        %Fsol=[0; log(1+(1-1/(L*delta*jMSY))^0.5)/log(2); log(1-(1-1/(L*delta*jMSY))^0.5)/log(2)]; %equation without E (faster)
        %Fsol=[0; log(1+(1-1/(L*delta*(MSYlevel(jMSY)-log(delta))))^0.5)/log(2); log(1-(1-1/(L*delta*(MSYlevel(jMSY)-log(delta))))^0.5)/log(2)]; %equation without E (faster)
        Fsol=[0; log(1+(1-1/(eL*(1/theta)*MSYlevel(jMSY)))^0.5)/log(2); log(1-(1-1/(eL*(1/theta)*MSYlevel(jMSY)))^0.5)/log(2)]; %equation without E (faster)
        Fsol(imag(Fsol)~=0)=-100;
        Fsol=sort(Fsol,'descend');
        Fsol(Fsol==-100)=NaN;
        Fsols(1:length(Fsol),i)=Fsol;
    end
    plot1=plot(Lrange,Fsols','Color','k','LineWidth',2);
    plot1(2).LineStyle='--';
    ylim([-6.5 3.5])
    title(['ln(MSY[kg])=' num2str(MSYlevel(jMSY))])
    zlabel 'log_2(F/F_M_S_Y)'
    xlabel 'cost/benefit (\gamma)'
    ylabel 'log_2(F/F_{MSY})'
    if jMSY==2
        legend 'F_{init}>F_{MSY}' 'F_{init}<F_{MSY}'
    end
end

%plot log2(F/Fmsy) vs. ln(MSY) with model predictions at I/V=1
%plot model fit at various MSYs
figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/3 scrsz(4)/2]);
colormap hsv
indicesBelowMSY=find(reshape(firstPt,1,1,[])<0);
indicesAboveMSY=find(reshape(firstPt,1,1,[])>=0);
%scatter3(L2,log(allMSY*1000),eF,[],reshape(log(allMSY*1000),[],1),'filled','jitter','on')
%scatter3(L2,log(allMSY*1000),eF,'k','filled','jitter','on')
scatter(log(allMSY(indicesAboveMSY)*1000),eF(indicesAboveMSY),'ok','filled','jitter','on','jitterAmount',0.02)
hold on
scatter(log(allMSY(indicesBelowMSY)*1000),eF(indicesBelowMSY),'ok','jitter','on','jitterAmount',0.02)
xlim([10 22.5])
ylim([-6.5 3.5])
grid on
hold on
eL=1;
lnMSYrange=[10:0.05:22];
for i=1:length(lnMSYrange)
    MSYlevel=lnMSYrange(i);
    Fsol=[0; log(1+(1-1/(eL*(1/theta)*MSYlevel))^0.5)/log(2); log(1-(1-1/(eL*(1/theta)*MSYlevel))^0.5)/log(2)];
    Fsol(imag(Fsol)~=0)=-100;
    Fsol=sort(Fsol,'descend');
    Fsol(Fsol==-100)=NaN;
    Fsols(1:length(Fsol),i)=Fsol;
end
plot1=plot(lnMSYrange,Fsols','Color','k','LineWidth',2);
plot1(2).LineStyle='--';
ylabel 'log_2(F/F_M_S_Y)'
xlabel 'ln(MSY[kg])'
legend 'F_{init}>F_{MSY}' 'F_{init}<F_{MSY}'
%title('fitted model')

%AIC
numPts=length(avgAllPt2);
LL_CASModel=-(numPts*log(2*pi)/2+numPts*log(std(eF))+LLweightedmeans/(2*var(eF))); %LLweightedmeans is actually sum of squares
SS_StatModel=(1-linstats(1))*STotmeansNull;
LL_StatModel=-(numPts*log(2*pi)/2+numPts*log(std(eF))+SS_StatModel/(2*var(eF)));
[AIC_CASModel BIC_CASModel]=aicbic(LL_CASModel,1,numPts)
[AIC_StatModel BIC_StatModel]=aicbic(LL_StatModel,4,numPts)
AICc_CASModel=AIC_CASModel+2*1*(1+1)/(numPts-1-1); %AIC+0.0186
AICc_StatModel=AIC_StatModel+2*4*(4+1)/(numPts-4-1); %AIC+0.189
diff_AIC=AIC_CASModel-AIC_StatModel
diff_AICc=AICc_CASModel-AICc_StatModel
diff_BIC=BIC_CASModel-BIC_StatModel

%plot bifurcation along L*ln(MSY) on x-axis
figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/3 scrsz(4)/2]);
colormap hsv
indicesBelowMSY=find(reshape(firstPt,1,1,[])<0);
indicesAboveMSY=find(reshape(firstPt,1,1,[])>=0);
%scatter3(L2,log(allMSY*1000),eF,[],reshape(log(allMSY*1000),[],1),'filled','jitter','on')
%scatter3(L2,log(allMSY*1000),eF,'k','filled','jitter','on')
scatter(L(indicesAboveMSY).*log(allMSY(indicesAboveMSY)*1000),eF(indicesAboveMSY),'ok','filled','jitter','on','jitterAmount',0.02)
hold on
scatter(L(indicesBelowMSY).*log(allMSY(indicesBelowMSY)*1000),eF(indicesBelowMSY),'ok','jitter','on','jitterAmount',0.02)
% xlim([10 22.5])
% ylim([-6.5 3.5])
grid on
hold on
%eL=1;
LlnMSYrange=[5:0.5:45];
Fsols=[];
for i=1:length(LlnMSYrange)
    LMSYlevel=LlnMSYrange(i);
    Fsol=[0; log(1+(1-1/((1/theta)*LMSYlevel))^0.5)/log(2); log(1-(1-1/((1/theta)*LMSYlevel))^0.5)/log(2)];
    Fsol(imag(Fsol)~=0)=-100;
    Fsol=sort(Fsol,'descend');
    Fsol(Fsol==-100)=NaN;
    Fsols(1:length(Fsol),i)=Fsol;
end
plot1=plot(LlnMSYrange,Fsols','Color','k','LineWidth',2);
plot1(2).LineStyle='--';
ylabel 'log_2(F/F_M_S_Y)'
xlabel '\gammaln(MSY[kg])'
legend 'F_{init}>F_{MSY}' 'F_{init}<F_{MSY}'
%title('fitted model')
