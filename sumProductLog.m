function uhat = sumProductLog(rx_cw, H, rho_dB, Niter)
% Log-domain sum product algorithm LDPC decoder
% Page 220, WILLIAM E. RYAN, SHU LIN, Channel Codes: Classical and Modern. 
% Inputs
%  rx_cw     : Received signal vector (column vector)
%  H         : LDPC matrix
%  N0        : Noise variance
%  Niter     : num of iterations
%
% Outputs
%  uhat      : Decoded vector (0/1)

[M, N] = size(H);
rho_lin = 10^(rho_dB/10);
N0 = 1/rho_lin;
% prior log-likelihood. Minus sign is used for 0/1 to -1/1 mapping
Lci = (-4*rx_cw./N0)';

% Init
LR = zeros(M, N);
phibetaij = zeros(M, N);
uhat = zeros(N,1);

% Asscociate the L(ci) matrix with non-zero elements of H
LQ = H.*repmat(Lci, M, 1);

% Get non-zero elements
[r, c] = find(H);

for n = 1:Niter

    % Get the sign and magnitude of L(Q)
    alphaij = sign(LQ);
    betaij = abs(LQ);

    for l = 1:length(r)
        % this comes from the (5.20)
        phibetaij(r(l), c(l)) = log((exp(betaij(r(l), c(l))) + 1)/...
                                   (exp(betaij(r(l), c(l))) - 1));
    end

    % loop over each check node
    for ii = 1:M

        % Find non-zeros in the row, i.e., find the connected variable nodes
        c1 = find(H(ii, :));

        % Get the summation of phi(betaij))
        for k = 1:length(c1)

            sumOfphibetaij = 0;
            prodOfalphaij = 1;

            % Summation of Pi(betaij)\c1[k]
            Pic1 = phibetaij(ii, c1);
            sumOfphibetaij = sum(Pic1) - Pic1(1,k);
%                sumOfPibetaij = sum(Pibetaij[i, c1]) - Pibetaij[i, c1[k]];

            % Avoid division by zero/very small number, get Pi(sum(Pi[betaij]))
            if sumOfphibetaij < 1e-20
                sumOfphibetaij = 1e-10;
            end
            phiSumOfPibetaij = log((exp(sumOfphibetaij)+1)/(exp(sumOfphibetaij) - 1));

            % Multiplication of alphaij\c1[k] (use '*' since alphaij are -1/1s)
            prodOfalphaij = prod(alphaij(ii, c1))*alphaij(ii, c1(k));

            % Update L[R]
            LR(ii, c1(k)) = prodOfalphaij*phiSumOfPibetaij;

        end % for k

    end % for ii

    % loop over each variable nodes
    for jj = 1:N

        % Find non-zeros in the column, i.e., find the connected check nodes
        r1 = find(H(:, jj));

        for k = 1:length(r1)

            % Update L[Q] by summation of L[rij]\r1[k]
            LQ(r1(k), jj) = Lci(jj) + sum(LR(r1, jj)) - LR(r1(k), jj);

        end % for k

        % Get L[Qi]
        LQi = Lci(jj) + sum(LR(r1, jj));

        % Decode L(Qi)
        if LQi < 0
            uhat(jj) = 1;
        else
            uhat(jj) = 0;
        end

    end % for jj
    
    % check if the decoding is already successful
    symdrome = mod(H*uhat, 2);
    if isempty(find(symdrome, 1))
        break
    end
        
end % for n