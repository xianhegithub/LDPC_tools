function out = J_sigma_inv(ei)
% this function is one of the steps in the EXIT calculation
% Ref. [1] Design of Low-Density Parity-Check Codes for Modulation and
% Detection -- Stephan ten Brink et al. Appendix

ei_star = 0.3646;
as1 = 1.09542;
bs1 = 0.214217;
cs1 = 2.33727;
as2 = 0.706692;
bs2 = 0.386013;
cs2 = -1.75017;

if ei >= 0 && ei <= ei_star
   
    out = sum([as1 bs1 cs1].* (ei.^([2 1 0.5])));
    
elseif ei > ei_star && ei < 1
    
    out = - as2 * log(bs2*(1 - ei)) - cs2 * ei;
    
else
    
    out = 10;
    fprintf('Numerical error in the inverse J_sigma function\n');
    
end