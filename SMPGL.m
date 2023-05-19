function [ x, cost ] = SMPGL( y, B, lam, rho, mu, Nit)

if nargin < 6
    Nit = 25;       % Default value
end
if nargin < 5   
    L = size(B, 2);
    mu = 2/L - 0.01;        % Default value
end
if nargin < 4
    rho = 0.0001;   % Default value
end
y = y(:);
L = size(B, 2);
cost = zeros(Nit , 1);
M = sum(B{1,1});



for i = 1 : L
    x(:, i) = y;
end
w = zeros(size(x));
v = w;
Tmp = w;
T = w;

for iter = 1:Nit
    r_tmp = [];
    for i = 1 : L
        r = sqrt(conv(abs(x(:, i)).^2, B{1, i}));
        r_tmp = [r_tmp, sum(r)];
        w(:, i) = 1 ./ (abs(x(:, i)) + eps);
        rr = conv(1./(r + eps), B{1, i} , 'valid');
        v(:, i) = 1 + lam * mu * rr;
    end      

    Tmp_sum = sum(x, 2) - y;
    for i = 1 : L
        Tmp(:, i) = (x(:, i) - mu * (Tmp_sum)) ./ (v(:, i) + eps);
        T(:, i) = mu * lam * rho * M ./ (v(:, i)) .* w(:, i);
        x(:, i) = softth(Tmp(:, i) , T(:, i));
    % enhance sparsity by reweighting
    end
    cost(iter) = 0.5 * sum(abs(y - sum(x, 2)).^2) + lam * sum(r_tmp + rho * M * sum(abs(x) ./ w, 1));
end

end