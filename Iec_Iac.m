function out = Iec_Iac(I_ac,  dc)
% this function calculate the EXIT curve for check nodes with degree dc
% Ref. [1] Channel Codes classical and modern -- William E. Ryan and Shu Lin 
% (9.43)

I_ec = zeros(1,length(I_ac));

for ii = 1:length(I_ac)
    
    tmp = J_sigma_inv(1 - I_ac(ii));
    J_arg = sqrt((dc - 1)*tmp^2);
    I_ec(ii) = 1 - J_sigma(J_arg);
    
end

out = I_ec;