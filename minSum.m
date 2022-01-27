function uhat = minSum(rx_cw, H, Niter)
% Simplified log-domain sum product algorithm LDPC decoder, known as the
% min-sum decoder.
% Page 226, Sec. 5.5.1, WILLIAM E. RYAN, SHU LIN, Channel Codes: Classical and Modern. 
% Inputs
%  rx_cw     : Received signal vector (column vector)
%  H         : LDPC matrix
%  Niter     : num of iterations
%
% Outputs
%  uhat      : Decoded vector (0/1)

[M, N] = size(H);

% Prior log-likelihood (simplified). Minus sign is used for 0/1 to -1/1 mapping
Lci = -rx_cw';

% Init
LR = zeros(M, N);
uhat = zeros(N,1);

% Asscociate the L(ci) matrix with non-zero elements of H
LQ = H.*repmat(Lci, M, 1);

for n = 1:Niter


    % Get the sign and magnitude of L(Q)
    alphaij = sign(LQ);
    betaij = abs(LQ);

    % loop over each check node
    for ii = 1:M

        % Find non-zeros in the row, i.e., find the connected variable nodes
        c1 = find(H(ii, :));

        % Get the minimum of betaij
        for k = 1:length(c1)

            % Minimum of betaij\c1[k]
            minOfbetaij = realmax();
            for l = 1:length(c1)
                if l ~= k
                    if betaij(ii, c1(l)) < minOfbetaij
                        minOfbetaij = betaij(ii, c1(l));
                    end
                end
            end % for l

            % Multiplication alphaij\c1(k) (use '*' since alphaij are -1/1s)
            prodOfalphaij = prod(alphaij(ii, c1))*alphaij(ii, c1(k));

            % Update L[R]
            LR(ii, c1(k)) = prodOfalphaij*minOfbetaij;

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

    end % for j

    % check if the decoding is already successful
    symdrome = mod(H*uhat, 2);
    if isempty(find(symdrome, 1))
        break
    end    
    
end % for n