% this script runs an example of constructing parity-check matrix of LDPC
% code using the progressive edge growth algorithm.

n = 4088; % number of variable nodes
m = 2044; % number of check nodes
% bit node degree polynomial
lambda = zeros(1,10);
lambda(2) = 0.451; lambda(3) = 0.3931; lambda(10) = 0.1558; 
% check node degree polynomial 
rho = zeros(1,8);
rho(8) = 0.7206; rho(9) = 0.2994;
code_rate = 1 - (sum(rho./(1:length(rho))))/(sum(lambda./(1:length(lambda))));

H = zeros(m, n);
v_deg = ceil(n * lambda);
if sum(v_deg) ~= n
    fprintf('v_deg has an error\n');
end

dv = [];
tmp = find(v_deg);
for jj = 1:length(tmp)
    dv = [dv ones(1, v_deg(tmp(jj)))*tmp(jj)];
end
dv = dv(randperm(n));

H = PEG_ACE(n, m, dv, 0);