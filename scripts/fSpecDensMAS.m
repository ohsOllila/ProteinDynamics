 function J = fSpecDensMAS(omega,omegaR,Itauc,Itaus)

J = 1/3*(fLorentzSpecDens(omega-omegaR,Itauc+Itaus)+fLorentzSpecDens(omega+omegaR,Itauc+Itaus)) + ...
    1/6*(fLorentzSpecDens(omega-2*omegaR,Itauc+Itaus)+fLorentzSpecDens(omega+2*omegaR,Itauc+Itaus));
