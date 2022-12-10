%% "Experiment_1.m" will generate Figs. 7-9 in the paper: 
%% Zheng Cao, Jisheng Dai, Weichao Xu, and Chunqi Chang, "Sparse Bayesian Learning Approach for Compound Bearing Fault Diagnosis"

clear;
close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%%  Paderborn University Dataset
 y=load('N15_M01_F10_KB27_10.mat');
 y=y.N15_M01_F10_KB27_10.Y(7).Data; 
 y=y';
 y=y(1:64000);
 Fs=64000;                                        %  the sampling frequency             
 Sig_N=y;
 N=length(Sig_N);
 t = (0 : N-1) / Fs;
 
 F1=76.25;                                        %   BPFO
 F2=123.3;                                        %   BPFI

y_h= abs(hilbert(Sig_N));                         %   the envelope of Sig_N
F = ([1:N]-1)*Fs/N;                               %   frequency domain

%%   plot original signal and its envelope spectrum
figure(6);
subplot(3,2,1)
plot(t,Sig_N,'black')
axis([0 1 -2 2])
xlabel('Time [s]')
ylabel('Amp.[m/s^2]')
title('a) Original signal')

subplot(3,2,2)
x_plot=[F1  F1];
y_plot=[0  0.03];
plot(x_plot,y_plot,'--g','linewidth',1);
hold on;
x_plot2=[F2  F2];
y_plot2=[0  0.03];
plot(x_plot2,y_plot2,'--r','linewidth',1);
hold on;
plot(F, abs(fft(y_h))/(N/2));
axis([0 800 0 0.03])
xlabel('Frequency [Hz]')
title('b) Envelope spectrum')

%% Our method
Fb=[F1, F2];
mu_SBL=Compound_fault_learning(Sig_N, Fs, Fb);
y_our1=abs(fft(abs(hilbert(mu_SBL(:,1))) -mean(abs(hilbert(mu_SBL(:,1)   ))) ))/(N/2);
y_our2=abs(fft(abs(hilbert(mu_SBL(:,2) )) -mean(abs(hilbert(mu_SBL(:,2) )))  ))/(N/2);

figure(7);
F_area= F(1:2001);
subplot(3,2,1)
plot(t,mu_SBL(:,1),'black')
axis([0 1 -0.2 0.2])
ylabel('Amp.[m/s^2]')
title('a) Our method: Impluse Outer')

subplot(3,2,2)
plot(F_area, y_our1(1:2001) );
axis([0 600 0 0.01])
title('b) Our method: Envelope Outer')

subplot(3,2,3)
plot(t,mu_SBL(:,2),'black')
axis([0 1 -0.2 0.2])
xlabel('Time [s]')
ylabel('Amp.[m/s^2]')
title('c) Our method: Impluse Inner')

subplot(3,2,4)
plot(F_area,  y_our2(1:2001))
axis([0 1000 0 0.01])
xlabel('Frequency [Hz]')
title('d) Our method: Envelope Inner')

%%   Set the parameters of SMPGL
% Fault1
K1 = 1;
N1 = 4;
N0 = round(Fs/F1) - N1;   
M = 4;
B1 = binaryblock( K1 , N0 , N1 , M );
% Fault2
K1 = 1;
N1 = 4;
N0 = round(Fs/F2) - N1;   
M = 4;
B2 = binaryblock( K1 , N0 , N1 , M );

%%   Perform the SMPGL algorithm
[C,L]=wavedec(Sig_N,5,'sym8');
c1=detcoef(C,L,1);
est_noise=median(abs(c1-median(c1)))/0.6745;
[lambda, rho] = Selection_Pars(est_noise);
B = {B1,B2};
Nit = 50;
mu = 2/3*0.9;
[x , cost1] = SMPGL(Sig_N , B, lambda , rho, mu, Nit);       
x_1 = x(:, 1);
x_2 = x(:, 2);
y_SMPGL1=abs(fft(abs(hilbert(x_1)) -mean(abs(hilbert(x_1)))  ))/(N/2);
y_SMPGL2=abs(fft(abs(hilbert(x_2)) -mean(abs(hilbert(x_2)))  ))/(N/2);

figure(8);
subplot(3,2,1)
plot(t,x_1,'black')
axis([0 1 -0.5 0.5])
ylabel('Amp.[m/s^2]')
title('a) SMPGL: Impluse Outer')

subplot(3,2,2)
plot(F_area, y_SMPGL1(1:2001) );
axis([0 600 0 0.015])
title('b) SMPGL: Envelope Outer')

subplot(3,2,3)
plot(t,x_2,'black')
axis([0 1 -0.5 0.5])
xlabel('Time [s]')
ylabel('Amp.[m/s^2]')
title('c) SMPGL: Impluse Inner')

subplot(3,2,4)
plot(F_area,  y_SMPGL2(1:2001))
axis([0 1000 0 0.015])
xlabel('Frequency [Hz]')
title('d) SMPGL: Envelope Inner')
