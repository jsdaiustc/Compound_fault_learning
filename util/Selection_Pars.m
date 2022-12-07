function [ lambda , rho ] = Selection_Pars(NoiseSigma)

load data_N1_4_M_4_Final.mat
lambda = coef(1) * NoiseSigma + coef(2);
if NoiseSigma > Sigma(N)
    rho = First;
else
    rho = coef_rho(1) * NoiseSigma.^2 + coef_rho(1) * NoiseSigma + coef_rho(3);
end

if rho > 0.01
    rho = 0.01;
end


end