clear
close all
% this script calculate the EXIT chart for irregular LDPC ensemble
% Ref. [1] Channel Codes classical and modern -- William E. Ryan and Shu Lin
% (9.44) (9.45)


% bit node degree polynomial
lambda = zeros(1,10);
lambda(2) = 0.267; lambda(3) = 0.176; lambda(4) = 0.127; lambda(10) = 0.43;
% check node degree polynomial
rho = zeros(1,8);
rho(5) = 0.113; rho(8) = 0.887;
code_rate = 1 - (sum(rho./(1:length(rho))))/(sum(lambda./(1:length(lambda))));

Ps = 1;
EbN0_dB = 0.55;
EbN0 = 10^(EbN0_dB/10);
sigma = sqrt(Ps/EbN0);
sigma_ch = sqrt(8 * code_rate * EbN0);


I_av = 0:0.05:1;
I_ev = zeros(1,length(I_av));

for ii = 1: length(lambda)

    I_ev = lambda(ii) * Iev_Iav(I_av, sigma_ch, ii) + I_ev;

end

plot(I_av, I_ev)

I_ac = 0:0.05:1;
I_ec = zeros(1,length(I_av));

for ii = 1: length(rho)

    I_ec = rho(ii) * Iec_Iac(I_ac, ii) + I_ec;

end

hold on
plot(I_ec, I_ac)
