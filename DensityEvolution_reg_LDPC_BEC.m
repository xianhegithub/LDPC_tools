% script that reproduces the result on page 43 of 
% Introducing Low-Density Parity-Check Codes -- Sarah J. Johnson
% open ~/Documents/trello/ldpctutorial.pdf

%era_prob = [0.3 0.42 0.43 0.5];
era_prob = [0.4293 0.4294]; 
iter_max = 300;
for ep = 1:length(era_prob)
wc = 3;
wr = 6;
iter = 7;
p0 = era_prob(ep);
rx_era_prob = zeros(1, iter_max);

for jj = 1:iter_max
    
    iter = jj;
    p = zeros(1,iter);
    
    for ii = 1:iter

        if ii == 1
            p(1) = p0 * (1 - (1 - p0)^(wr - 1))^(wc - 1);
        else
            p(ii) = p0 * (1 - (1 - p(ii-1))^(wr - 1))^(wc - 1);
        end

    end
    rx_era_prob(jj) = p(end);

end

plot(1:iter_max, rx_era_prob)
hold on
end
%legend('p0=0.3','p0=0.42','p0=0.43','p0=0.5')
legend('p0=0.4293','p0=0.4294')