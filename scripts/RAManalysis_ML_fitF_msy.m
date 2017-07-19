function LL = RAManalysis_ML_fitF_msy(avgAllPt,firstPt,params)
digits(5); %decrease vpa precision
syms u(F) F

theta=params(1);
E=params(2);

%V=avgAllPt(:,3,:); %mean landing value/landings over years
eF=avgAllPt(:,1,:); %mean harvest rate over years
L=avgAllPt(:,2,:); %mean (cost+subsidy)/landing value over years
%MSY=log10(avgAllPt(:,5,:)+1);
MSY=log(avgAllPt(:,5,:)*E);
%MSY=avgAllPt(:,5,:)*E;
V=1; %set V to be constant

%evaluate fit of fishery mean harvest rate data to modelled equilibria
LLmeans=[];
for i=1:length(L)
    %for full model:
    %     u(F)=vpa(delta*MSY(i)*(1-F).*V.*(1/(delta*MSY(i)*F.*(1-F/2))-L(i))+E); %equation with E
    %     Fsol=vpasolve(eval(u(F)),F);
    %     Fsol(imag(Fsol)~=0)=-100;
    %     Fsol(Fsol<0)=-100;
    %     Fsol=sort(Fsol,'descend');
    %     Fsol(Fsol==-100)=NaN;
    %     Fsol=log2(Fsol);
    
    %for model without imperfection:
    %Fsol=[0; log(1+(1-1/(L(i)*delta*MSY(i)))^0.5)/log(2); log(1-(1-1/(L(i)*delta*MSY(i)))^0.5)/log(2)]; %equation without E (faster)
    Fsol=[0; log(1+(1-1/(L(i)*(1/theta)*MSY(i)))^0.5)/log(2); log(1-(1-1/(L(i)*(1/theta)*MSY(i)))^0.5)/log(2)]; %equation without E (faster)
    %Fsol=[0; log(1+(1-1/(L(i)*log(E*MSY(i)-delta)+delta*E*MSY(i)))^0.5)/log(2); log(1+(1-1/(L(i)*log(E*MSY(i)-delta)+delta*E*MSY(i)))^0.5)/log(2)]; %equation without E (faster)
    Fsol(imag(Fsol)~=0)=-100;
    Fsol=sort(Fsol,'descend');
    Fsol(Fsol==-100)=NaN;
    
    
    %if eF(i)<Fsol(2) %mean harvest rate of fishery is less than the unstable equilibrium
    if firstPt(i)<Fsol(2) %first known harvest rate of fishery is less than the unstable equilibrium, thus closer to second lower stable equilibrium (and there is one)
        %      if ~isnan(Fsol(3))
        LLmeans(i)=(eF(i)-Fsol(3))^2;
        %      else
    else %datapoint is closer to first higher stable equilibrium (or the model only predicts single equilibrium)
        LLmeans(i)=(eF(i)-Fsol(1))^2;
    end
    if isnan(LLmeans(i))
        LLmeans(i)=Inf;
    end
end
LL=sum(LLmeans);