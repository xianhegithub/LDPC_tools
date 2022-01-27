function out = J_sigma(sigma)
% this function is one of the steps in the EXIT calculation
% Ref. [1] Design of Low-Density Parity-Check Codes for Modulation and
% Detection -- Stephan ten Brink et al. Appendix


sigma_star = 1.6363;
aj1 = -0.0421061;
bj1 = 0.209252;
cj1 = -0.00640081;
aj2 = 0.00181491;
bj2 = -0.142675;
cj2 = -0.0822054;
dj2 = 0.0549608;

if sigma >= 0 && sigma <= sigma_star
    
    out = sum([aj1 bj1 cj1].* (sigma.^([3 2 1])));
    
elseif sigma > sigma_star
    
    out = 1 - exp(sum([aj2 bj2 cj2 dj2].* (sigma.^([3 2 1 0]))));
    
else
    
    out = 1;
    
end
   