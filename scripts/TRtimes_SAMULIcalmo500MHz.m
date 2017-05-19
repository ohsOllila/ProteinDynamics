% THE CORRELATION FUNCTION TO BE ANALYZED SHOULD BE IN corrF
%
%   load rotaCF_10.dat
%   corrF=rotaCF_10;
% %function rates=FT(corrF)
%
%

%NCtimes=3;
%Ctimes=(-10:0.1:NCtimes);
%Ctimes=0.1*10.^Ctimes;
Ctimes=(0.01:0.1:50);
Ctimes=Ctimes*0.01*10^(-9);

% SET ORDER PARAMETER TO ZERO FOR LIQUID SAMPLE
%carbon=1;
%ordPcomp(1)=0;
%Ccoeffs=(1:100);
%load ordPcomp_simTMP.dat
%Ctimes=Ctimes*0.001*10^(-9);

% NORMALIZED CORRELATION FUNCTION
NcorrF=(corrF-ordP^2)/(1-ordP^2);

%%%%%UNFORCED PLATEAU: CORRELATION FUNCTION NOT FORCED TO DECAY TO PLATEAU
for i = 1:(length(NcorrF))
    for j = 1:1:length(Ctimes)
        Cexp(i,j)=exp(-NcorrF(i,1)*0.01*10^(-12)/(Ctimes(j)));
    end
    %Cexp(i+1,length(Ctimes)+1)=1;
end
Coeffs=lsqnonneg(Cexp,NcorrF(:,2));



%%%%FORCED PLATEAU: CORRELATION FUNCTION FORCED TO REACH THE S^2 (ZERO IN LIQUID)
for i = 0:1:length(corrF)-1
    for j = 1:1:length(Ctimes)
        %S=0.006;
        Cexp2(i+1,j)=(1-(ordP)^2)*exp(-i*10^(-12)/(Ctimes(j)))+(ordP)^2;
      %         Cexp(i+1,j)=exp(-i/(Ctimes(j)));

    end
    %Cexp(i+1,length(Ctimes))=(ordPcomp(carbon))^2;
end
Coeffs2=lsqnonneg(Cexp2,corrF(:,2));
% 
%  
%    figure(1),clf
%  plot(corrF(:,1),corrF(:,2))
% plot(corrF(:,1),(Cexp(:,1:length(Ctimes))*Coeffs(1:length(Ctimes)))+Coeffs(length(Ctimes)+1),'b')
% %   plot(corrF(:,1),(Cexp(:,1:length(Ctimes))*Coeffs2(1:length(Ctimes)))+Coeffs2(length(Ctimes)+1),'r')
% return
  % % plot(corrF(:,1),(Cexp(:,1:length(Ctimes))*Coeffs(1:length(Ctimes))),corrF(:,1),Cexp*Coeffs,'r');
%  return

%Ctimes=Ctimes*0.001*10^(-9);
% 
%
%MAKE FITTED FUNCTIONS FOR PLOTTING
fitPL(1:2*length(corrF),2)=0;
fitPL(1:length(corrF),1)=corrF(:,1);
fitPL(length(corrF)+1:2*length(corrF),1)=corrF(:,1)+corrF(length(corrF),1);

fitPL2(1:length(corrF),2)=0;
fitPL2(:,1)=corrF(:,1);

% %%%%UNFORCED PLATEAU
 for i = 1:1:2*length(NcorrF)
     for j = 1:1:length(Ctimes)
         %fitPL(i,2)=fitPL(i,2)+Coeffs(j)*((1-ordP^2)*Cexp(i,j)+ordP^2);  %Coeffs(j)*exp(-fitPL(i,1)/(Ctimes(j))); 
     end
 end

%%%%%%FORCED PLATEAU0
for i = 1:1:length(corrF)
    for j = 1:1:length(Ctimes)
        fitPL2(i,2)=fitPL2(i,2)+Coeffs2(j)*Cexp2(i,j);
    end
end

% %fitPL
%subplot(1,2,1);
%figure(1),clf
%plot(fitPL(:,1),fitPL(:,2),'b')
%   hold on
%  plot(corrF(:,1),corrF(:,2),'r')
%,corrF(:,1),Coeffs(length(Ctimes)+1)*ones(length(corrF),1)
% 
% %subplot(1,2,2);
%  plot(fitPL2(:,1),fitPL2(:,2),'b')
%    hold on
%   plot(corrF(:,1),corrF(:,2),'r',corrF(:,1),(ordPcomp_sim(carbon))^2*ones(length(corrF),1))

 %  %fitPL2
%figure(2),clf
% subplot(1,2,1);
% plot(fitPL2(:,1),fitPL2(:,2),'b',fitPL(:,1),(ordPcomp(carbon))^2*ones(length(fitPL),1),'b--')
 %  hold on
 % plot(corrF(:,1),corrF(:,2),'r')
 % hold on

%SET TIME SCALE FOR X-AXIS
%Ctimes=Ctimes*0.001*10^(-9);
 
%EFFECTIVE CORRELATION TIME
tau_eff=0;
for i=1:length(Ctimes)
    tau_eff=tau_eff+(Coeffs(i)*Ctimes(i));
end

area1_length=length(NcorrF(:,2));
area1=(NcorrF(area1_length,2)*NcorrF(area1_length,1)*10^-12);
area=(sum(NcorrF(:,2))*(NcorrF(2,1)-NcorrF(1,1))*10^-12);%-area1;
area2=(sum(fitPL(:,2))*(NcorrF(2,1)-NcorrF(1,1))*10^-12);

% tau_eff2=0;
% 
% for i=1:length(Ctimes)
%     tau_eff2=tau_eff2+(Coeffs(i)*Ctimes(i));
% end


% 
% fitPLft=numFTINT(fitPL);
% corrFft=numFTINT(corrF);
% subplot(1,2,2);
% plot(log10(fitPLft(:,1)/(0.001*10^(-9))),log10(0.001*10^(-9)*fitPLft(:,2)));
% plot((fitPLft(:,1)/(0.001*10^(-9))),2*3.14*(0.001*10^(-9)*fitPLft(:,2)));
% hold on
% plot(log10(corrFft(:,1)/(0.001*10^(-9))),log10(0.001*10^(-9)*corrFft(:,2)),'r')
% plot((corrFft(:,1)/(0.001*10^(-9))),(0.001*10^(-9)*corrF(:,2)),'r')
% return


%CALCULATE SPECTRAL DENSITY
JhMn=0;
Jn=0;
JhPn=0;
J0=0;
Jh=0;
Jw1=0;

J12=0;

J(1:length(Ctimes),1:2)=0;
J(1,1)=Coeffs(length(Coeffs))*0.001*10^(-9);

for j= 1:1:8000
    SpectDfit(j,2)=0;
end

for j= 1:1:8000
    w=(j-1)*10*3.14*10^9/2000; 
    SpectDfit(j,1)=w;
    for i = 1:1:length(Ctimes)
        SpectDfit(j,2)=SpectDfit(j,2)+Coeffs(i)*fLorentzSpecDensSAMULI(w,Ctimes(i));   
    end
end

%subplot(1,2,2);
%plot(log10(SpectDfit(:,1)),log10(SpectDfit(:,2)),'b')
%plot((SpectDfit(:,1)),(SpectDfit(:,2)),'g')


SpectDfit(1,2);
SpectDfit(201,2);
SpectDfit(401,2);

%omegaR=2*3.14*5*10^(3);
omegaR=0;
%Itaus=1/(2*3.14*5*10^(4));

%Itaus=6.0*10^(-9);
Itaus=0;
%omegaR=0;

%JhMn=Coeffs2(length(Coeffs2))*fSpecDensMAS(0,omegaR,0,Itaus);
%JhMn=ordPcomp_simHGgly(carbon)^2*fSpecDensMAS(0,omegaR,0,Itaus);

%wh=500*10^(6);
%wc=10000;
%wcw=0.25*wh;
gammaH=267.513*10^6;
gammaC=67.262*10^6;
gammaN=-27.166*10^6;
wh=500*10^6;
wc=wh*gammaC/gammaH;
%wc=125.76*10^6;
wn=wh*gammaN/gammaH;
%wh=wc/0.25;
w1=50000;


 for i = 1:1:length(Ctimes)
        Ctimes(i)=100*Ctimes(i);
        w=wh-wn;
        %JhMn=JhMn+2*Coeffs(i)*Ctimes(i)/(1+w*w*Ctimes(i)*Ctimes(i));
        JhMn=(JhMn+Coeffs2(i)*fSpecDensMAS(w,omegaR,Ctimes(i),Itaus));
        w=wn;
        %Jn=Jn+2*Coeffs(i)*Ctimes(i)/(1+w*w*Ctimes(i)*Ctimes(i));
        Jn=(Jn+Coeffs2(i)*fSpecDensMAS(w,omegaR,Ctimes(i),Itaus));
        w=wh;
        Jh=(Jh+Coeffs2(i)*fSpecDensMAS(w,omegaR,Ctimes(i),Itaus));
        w=wn+wh;
        %JhPn=JhPn+2*Coeffs(i)*Ctimes(i)/(1+w*w*Ctimes(i)*Ctimes(i));
        JhPn=(JhPn+Coeffs2(i)*fSpecDensMAS(w,omegaR,Ctimes(i),Itaus));
        w=0;
        J0=(J0+Coeffs2(i)*fSpecDensMAS(w,omegaR,Ctimes(i),Itaus));
        w=w1;
        Jw1=(Jw1+Coeffs2(i)*fSpecDensMAS(w,omegaR,Ctimes(i),Itaus));
 end
 
mu=4*3.14*10^(-7);
h_planck=1.055*10^(-34);
r=0.109*10^(-9);
rN=0.101*10^(-9);
d=(mu*gammaN*gammaH*h_planck)/(4*pi*rN^3);
const=(((2*3.14*mu*gammaC*gammaH*h_planck)/(8*3.14*r^3))^2)/20;
b=1.3*10^(-3);
R1=(d^2/(4*5))*(JhMn+3*Jn+6*JhPn)+Jn*(2*pi*wn*160*10^(-6))^2/15;
R1_agn=2*(gammaC^2)*(b^2)*Jn;
R2=0.5*(d^2/(4*5))*(JhMn+3*Jn+6*JhPn+4*J0+6*Jh)+(4*J0+3*Jn)*(2*pi*wn*160*10^(-6))^2/(15*6);
NOE=1+(d^2/(4*5))*(6*JhPn-JhMn)*gammaH/(gammaN*R1);
%rates

%end