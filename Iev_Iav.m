function out = Iev_Iav(I_av, sigma_ch, dv)
% this function calculate the EXIT curve for variable nodes with degree dv
% Ref. [1] Channel Codes classical and modern -- William E. Ryan and Shu Lin 
% (9.42)

I_ev = zeros(1,length(I_av));
for ii = 1:length(I_av)
    
    tmp = J_sigma_inv(I_av(ii));
    J_arg = sqrt((dv - 1)*tmp^2+sigma_ch^2);
    I_ev(ii) = J_sigma(J_arg);
    
end

out = I_ev;