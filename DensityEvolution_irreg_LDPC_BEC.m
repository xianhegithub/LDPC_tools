% script that reproduces the result on page 43 of 
% Introducing Low-Density Parity-Check Codes -- Sarah J. Johnson
% open ~/Documents/trello/ldpctutorial.pdf
clear all
era_prob = [0.4 0.425 0.45 0.475 0.5];
%era_prob = [0.3 0.4 0.5 0.6];
%era_prob = [0.455 0.46 0.465 0.47 0.475];
iter_max = 50;
% bit node degree polynomial
lambda = zeros(1,20);
lambda(2) = 0.1; lambda(3) = 0.4; lambda(20) = 0.5; 
% check node degree polynomial
rho = zeros(1,9);
rho(8) = 0.5; rho(9) = 0.5;
code_rate = 1 - (sum(rho./(1:length(rho))))/(sum(lambda./(1:length(lambda))));


for ep = 1:length(era_prob)

p0 = era_prob(ep);
rx_era_prob = zeros(1, iter_max);

for jj = 1:iter_max
    
    iter = jj;
    p = zeros(1,iter);
    
    for ii = 1:iter

        if ii == 1
            p(1) = p0 * sum(lambda.*(1 - sum(rho.* (1 - p0).^(1:length(rho)))).^(1:length(lambda)));
        else
            p(ii) = p0 * sum(lambda.*(1 - sum(rho.* (1 - p(ii-1)).^(1:length(rho)))).^(1:length(lambda)));
        end

    end
    rx_era_prob(jj) = p(end);

end

plot(1:iter_max, rx_era_prob)
hold on
end
ylim([0 0.6])