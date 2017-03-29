clear all

load ~/Dropbox/Calmodulin/Data/residues.dat
%OPs=load('~/Dropbox/PsTonB/Data/OPfromPLATEAU.dat');

figure(2),clf
for ac = 1:1:144
    filename=['~/Dropbox/Calmodulin/Data/NHcorrelationFunctions/withAMBER_400-1000ns/scaledrotation/NHrotaCF_' num2str(ac) '.xvg']  
    corrFtmp=load(filename);
    corrF=[corrFtmp(1:60,1),corrFtmp(1:60,2)];
    %ordP=sqrt(OPs(ac,2));
    ordP=0;
    TRtimes_SAMULIcalmo
    R_1(ac+1)=R1;
    R_2(ac+1)=R2;
    J_0(ac+1)=J0;
    NOEs(ac+1)=NOE;
    p_sim2(:,ac+1)=Coeffs2;
    p_sim(:,ac+1)=Coeffs;
    CoeffsSAVED(:,ac+1)=Coeffs2;
    plot(10^(-3)*fitPL2(:,1),fitPL2(:,2),'b',fitPL(:,1),(ordP)^2*ones(length(fitPL),1),'b--')
    %plot(fitPL(:,1),fitPL(:,2),'b')
    %plot(10^(-3)*fitPL(:,1),(1-ordPcomp(carbon)^2)*fitPL(:,2)+ordPcomp(carbon)^2,'b',10^(-3)*fitPL(:,1),(ordPcomp(carbon))^2*ones(length(fitPL),1),'b--')
    axis([0 10^(-3)*fitPL(length(fitPL(:,1)),1)*2/3 0 1])
    xlabel('t (ns)')
    title('\beta')
    hold on
    plot(10^(-3)*corrF(:,1),corrF(:,2),'r')
    %plot(10^(-3)*NcorrF(:,1),(1-ordP^2)*NcorrF(:,2)+ordP^2,'r')
    hold on
end

figure(3),clf
plot(residues,1./R_1,'bo')
hold on
plot(residues,1./R_1)
load  ~/Dropbox/Calmodulin/Data/experimentalRELAXATIONdata/T1experimentalDATA.dat
plot(T1experimentalDATA(:,1),T1experimentalDATA(:,2))
xlabel('N number');
ylabel('T_1 / s^{-1}');
title('T_1');
T_1=[residues,1./R_1'];
T1file=fopen('~/Dropbox/Calmodulin/Data/T1fromSmulationsAMBERscaledrotation.dat','w');
fprintf(T1file,'%6.2f %12.8f\n',T_1');
fclose(T1file);
figure(4),clf
plot(residues,1./R_2,'bo')
hold on
plot(residues,1./R_2)
 %p=legend('harris','agnieszka',2);
load ~/Dropbox/Calmodulin/Data/experimentalRELAXATIONdata/T2experimentalDATA.dat
plot(T2experimentalDATA(:,1),T2experimentalDATA(:,2))
xlabel('N number');
ylabel('T_2 / s^{-1}');
title('T_2');
T_2=[residues,1./R_2'];
T2file=fopen('~/Dropbox/Calmodulin/Data/T2fromSmulationsAMBERscaledrotation.dat','w');
fprintf(T2file,'%6.2f %12.8f\n',T_2');
fclose(T2file);
figure(5),clf
NOEresults=[residues,NOEs']
plot(NOEresults(:,1),NOEresults(:,2))
hold on
load ~/Dropbox/Calmodulin/Data/experimentalRELAXATIONdata/15N-NOEexperimentalDATA.dat
plot(X15N_NOEexperimentalDATA(:,1),X15N_NOEexperimentalDATA(:,2))
NOEfile=fopen('~/Dropbox/Calmodulin/Data/NOEfromSmulationsAMBERscaledrotation.dat','w');
fprintf(NOEfile,'%6.2f %12.8f\n',NOEresults');
fclose(NOEfile);
figure(6),clf
plot(Ctimes'*10^9,CoeffsSAVED(:,:))
%plot(194:1:280,J_0*10^9,'bo')
%p=legend('harris','agnieszka',2);
%xlabel('N number');
%ylabel('J(0) ns/rad');
%title('J0');
