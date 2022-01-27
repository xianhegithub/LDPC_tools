clear
close all
% this script calculate the EXIT chart for regular LDPC ensemble 
% Ref. [1] Channel Codes classical and modern -- William E. Ryan and Shu Lin 
% (9.42) (9.43)


dv = 3;
dc = 6;
code_rate = dv/dc;
ll = 0;
Ps = 1;
EbN0_dB = 1.1;
EbN0 = 10^(EbN0_dB/10);
sigma = sqrt(Ps/EbN0);
sigma_ch = sqrt(8 * code_rate * EbN0);


I_av = 0:0.05:1;
I_ev = Iev_Iav(I_av, sigma_ch, dv);


plot(I_av, I_ev)

I_ac = 0:0.05:1;
I_ec = Iec_Iac(I_ac, dc);

hold on
plot(I_ec, I_ac)