clear
close all
% this script does density evolution for regular LDPC code ensemble in AWGN
% channel using Gaussian approximation.
% Ref. [1] Channel Codes classical and modern -- William E. Ryan and Shu Lin 
% (9.32) (9.36) and Example 9.6


iter_max = 10000;
Pe_th = 1e-6;

dv = 3;
dc = 6;
ll = 0;
sigma = 0.6;
y = -20:0.1:120;


while ll < iter_max
    
    
    
    if Pe < Pe_th
        fprintf('Pe threshold satified\n');
        break;
    end
    
    ll = ll + 1;
    
    
end