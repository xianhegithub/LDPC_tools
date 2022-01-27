clear
sigma = 0:0.1:20;
J = zeros(1, length(sigma));
mubar = zeros(1, length(sigma));

for ii = 1:length(sigma)
    J(ii) = J_sigma(sigma(ii));
end

for ii = 1:length(sigma)
   sigmabar(ii) = J_sigma_inv(J(ii));
end


plot(sigma, J)
hold on
plot(sigmabar, J)
