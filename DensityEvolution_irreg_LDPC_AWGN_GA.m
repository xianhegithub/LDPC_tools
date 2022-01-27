clear
close all
% this script does density evolution for regular LDPC code ensemble in AWGN
% channel using Gaussian approximation.
% Ref. [1] Channel Codes classical and modern -- William E. Ryan and Shu Lin 
% (9.37) (9.38)



iter_max = 10000;
Pe_th = 1e-6;

% bit node degree polynomial
lambda = zeros(1,20);
lambda(2) = 0.23403; lambda(3) = 0.21242; lambda(6) = 0.1469; lambda(7) = 0.10284; lambda(20) = 0.30381;
% check node degree polynomial 
rho = zeros(1,9);
rho(8) = 0.71875; rho(9) = 0.28125;
code_rate = 1 - (sum(rho./(1:length(rho))))/(sum(lambda./(1:length(lambda))));

ll = 0;
sigma = 0.947;
y = -20:0.1:120;

mu_0 = 2/sigma^2 * ones(1, length(lambda));
mu_c = 0;
mu_c_delta = zeros(1, length(rho));

while ll < iter_max
    
    mu_v = mu_0 + ((1:length(lambda)) - 1) * mu_c;
    
        gmpdf = zeros(1,length(y)); 
    for ii = 1:length(mu_v)
        tmp = normpdf(y, mu_v(ii), sqrt(2*mu_v(ii)));
        gmpdf = gmpdf + lambda(ii)*tmp;
    end
    Pe = (y(2)-y(1))*sum(gmpdf(1:201));
    
    if Pe < Pe_th
        fprintf('Pe threshold satified\n');
        break;
    end
    
    ll = ll + 1;
    tmpv = zeros(1, length(mu_v));
    for ii = 1:length(tmpv)        
        tmpv(ii) = phi_mu(mu_v(ii));
    end
    
    for ii = 1:length(rho)
        tmpc = 1 - (1 - sum(lambda.* tmpv)) ^ (ii - 1);
        mu_c_delta(ii) = phi_mu_inv(tmpc);
    end
    
    mu_c = sum(rho.*mu_c_delta);
    
end