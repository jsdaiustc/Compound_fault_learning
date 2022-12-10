%% "Simulation_3.m" will generate Fig. 5 in the paper: 
%% Zheng Cao, Jisheng Dai, Weichao Xu, and Chunqi Chang, "Sparse Bayesian Learning Approach for Compound Bearing Fault Diagnosis"

clear 
close all
addpath(genpath(fileparts(mfilename('fullpath'))));
rand('seed',25)
randn('seed',25)

%% Creat the simulation signal
Fs = 12800;                                 %% the sampling frequency
N = 12800;                                  %% the length of compound signal

% fault_1
F1 = 37;                                    %% fault characteristic frequency 1
Sig1 = CreatSimulation(N , F1 , Fs);        %% creat impulse signal 1

% fault_2
F2=49;                                      %% fault characteristic frequency 2, this fault does not truly occur.

% fault_3
F3 = 57;                                    %% fault characteristic frequency 3
Sig3 = CreatSimulation(N , F3 , Fs);        %% creat impulse signal 3

Sig1=circshift(Sig1, randperm(1000,1) );
Sig3=circshift(Sig3, randperm(1000,1) );


%%  generate the compound fault
t = (0 : N-1) / Fs;
Sigma =.8;
True_signal=Sig1 + Sig3;
Sig_N = True_signal + Sigma * randn(N ,1);


%%  Perform our method
Fb=[F1, F2,F3];
mu_SBL=Compound_fault_learning(Sig_N-mean(Sig_N), Fs, Fb);


%% Set the parameters for SNPGL
% Fault1
K1 = 1;
N1 = 4;
N0 = round(Fs/F1) - N1;   % Fault1
M = 4;
B1 = binaryblock( K1 , N0 , N1 , M );
% Fault2
K1 = 1;
N1 = 4;
N0 = round(Fs/F2) - N1;   % Fault2
M = 4;
B2 = binaryblock( K1 , N0 , N1 , M );
% Fault3
K1 = 1;
N1 = 4;
N0 = round(Fs/F3) - N1;   % Fault3
M = 4;
B3 = binaryblock( K1 , N0 , N1 , M );


%% Perform the SMPGL
[lambda, rho] = Selection_Pars(Sigma);
B = {B1, B2, B3};
Nit = 50;
mu = 2/3*0.9;
[x , cost1] = SMPGL(Sig_N , B, lambda , rho, mu, Nit);       
x_1 = x(:, 1);
x_2 = x(:, 2);
x_3 = x(:, 3);

%% Plot the results
figure(1);
subplot(3,2,1)
plot(t,mu_SBL(:,1),'black')
axis([0 1 -1 1])
ylabel('Amp.[m/s^2]')
title('a) Our method (P1=1/37s)')

subplot(3,2,2)
plot(t, x_1, 'black')
axis([0 1 -1 1])
title('b) SMPGL (P1=1/37s)')

subplot(3,2,3)
plot(0,0,'black')
axis([0 1 -1 1])
ylabel('Amp.[m/s^2]')
title('c) Our method (P2=1/49s)')

subplot(3,2,4)
plot(t, x_2, 'black')
axis([0 1 -1 1])
title('d) SMPGL (P2=1/49s)')

subplot(3,2,5)
plot(t,mu_SBL(:,2), 'black')
axis([0 1 -1 1])
xlabel('Time [s]')
ylabel('Amp.[m/s^2]')
title('e) Our method (P3=1/57s)')

subplot(3,2,6)
plot(t, x_3, 'black')
axis([0 1 -1 1])
xlabel('Time [s]')
title('f) SMPGL (P2=1/57s)')
