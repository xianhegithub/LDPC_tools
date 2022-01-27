function out = phi_mu(mu)
% the function is used in density evolution for LDPC code ensemble using
% Gaussian approximation
% Ref. [1] Channel Codes classical and modern -- William E. Ryan and Shu Lin 
% Example 9.6


if (mu > 0 && mu < 10)
    
    out = exp(-0.4527*mu^0.86 + 0.0218);
   
elseif mu == 0
    
    out = 1;
    
else
    
    out = sqrt(pi/mu)*exp(-mu/4)*(1 - 11.5/7/mu);   % In reference it be 10/7/mu, but to make the function monotonic
    
end
