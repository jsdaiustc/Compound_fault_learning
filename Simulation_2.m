%% "Simulation_2.m" provides an example for Simulation 2 in the paper: 

clear 
close all
addpath(genpath(fileparts(mfilename('fullpath'))));

rand('seed',5)
randn('seed',5)
%% Creat the simulation signal
Fs = 12800;                                          %% the sampling frequency
N = 12800;                                           %% the length of compound signal
F=2000;

% fault_1
F1 = 37;                                             %% fault characteristic frequency 1
Sig1 = Generate_Simulation_noseed(Fs,F,N,F1);        %% creat impulse signal 1


% fault_2
F2=49;                                               %% fault characteristic frequency 2, this fault does not truly occur.
Sig2 = Generate_Simulation_noseed(Fs,F,N,F2);        %% creat impulse signal 3


% fault_3
F3 = 57;                                             %% fault characteristic frequency 3
Sig3 = Generate_Simulation_noseed(Fs,F,N,F3);        %% creat impulse signal 3

Sig1=circshift(Sig1, randperm(1000,1) );
Sig2=circshift(Sig2, randperm(1000,1) );
Sig3=circshift(Sig3, randperm(1000,1) );



%%  generate the compound fault
t = (0 : N-1) / Fs;
Sigma =0.8;
True_signal=Sig1+ Sig3;
Sig_N = True_signal + Sigma * randn(N ,1);

F_area = ([1:N]-1)*Fs/N; 
y_Sig=abs(fft(abs(hilbert(Sig_N)) -mean(abs(hilbert(Sig_N   ))) ))/(N/2);


%%  Perform our method
Fb=[F1,F2,F3];
[mu_SBL]=Compound_fault_learning(Sig_N-mean(Sig_N), Fs, Fb);



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
mse_SMPGL=(norm(True_signal - sum(x,2),'fro')/norm(True_signal,'fro') )^2;
x_1 = x(:, 1);
x_2 = x(:, 2);
x_3 = x(:, 3);


%% Plot the results
figure(5);
subplot(4,2,1)
plot(t,Sig_N,'black')
axis([0 1 -2.5 2.5])
xlabel('Time [s]')
ylabel('Amp.[m/s^2]')
title('a) ')

subplot(4,2,2)
x_plot=[F1  F1];
y_plot=[0  0.12];
plot(x_plot,y_plot,'--g','linewidth',1);
hold on;
x_plot2=[F3  F3];
y_plot2=[0  0.12];
plot(x_plot2,y_plot2,'--b','linewidth',1);
hold on;
plot(F_area(1:400), y_Sig(1:400), 'blue')
axis([0 400 0 0.12])
xlabel('Frequency [Hz]')
title('b) ')

subplot(4,2,3)
plot(t,mu_SBL(:,1),'black')
axis([0 1 -2.5 2.5])
xlabel('Time [s]')
ylabel('Amp.[m/s^2]')
title('c) ')

subplot(4,2,4)
plot(t, x_1, 'black')
axis([0 1 -2.5 2.5])
xlabel('Time [s]')
title('d) ')

subplot(4,2,5)
plot(0,0, 'black')
axis([0 1 -2.5 2.5])
xlabel('Time [s]')
ylabel('Amp.[m/s^2]')
title('e) ')

subplot(4,2,6)
plot(t, x_2, 'black')
axis([0 1 -2.5 2.5])
xlabel('Time [s]')
title('f)')

subplot(4,2,7)
plot(t,mu_SBL(:,2), 'black')
axis([0 1 -2.5 2.5])
xlabel('Time [s]')
ylabel('Amp.[m/s^2]')
title('g) ')

subplot(4,2,8)
plot(t, x_3, 'black')
axis([0 1 -2.5 2.5])
xlabel('Time [s]')
title('h) ')





