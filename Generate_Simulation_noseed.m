function [ Simulation] = Generate_Simulation_noseed(Fs,F,N,Fn)

Simulation = zeros(N , 1);
K = floor(Fs / Fn);
Num = floor(N / K);
slip = round(2000*(2*rand(Num,1)-1)*0.001); 
Interval = 30;
n = (0 : Interval-1) / Fs;
n = n';
for i = 1 : Num
    U = randi(3);
    Transient = zeros(Interval , 1); 
    for j = 1 : U
        A = randn(1) ;
        w = sqrt(10)*randn(1) + F;
        beta =sqrt(1) * randn(1);
        Transient = Transient + A*sin(2*pi*w*n + beta) .* exp(-900*n);
    end 
    if i == 1
        Simulation((i-1)*K+1:(i-1)*K+Interval) = Transient;
    else
        Simulation((i-1)*K+1+slip(i):(i-1)*K+slip(i)+Interval) = Transient;
    end
end

end

