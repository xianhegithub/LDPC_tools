function uhat = sumProduct(rx_cw, H, rho_dB, Niter)
% Probability-domain sum product algorithm LDPC decoder
% Page 250, WILLIAM E. RYAN, SHU LIN, Channel Codes: Classical and Modern.  
% Inputs
%   zdata     : Received signal vector (column vector)
%   H         : LDPC matrix
%   N0        : Noise variance
%   Niter     : num of iterations
%
% Outputs
%   uhat      : Decoded vector (0/1)
[M, N] = size(H);
rho_lin = 10^(rho_dB/10);
N0 = 1/rho_lin;

% priors
P1 = ones(size(rx_cw))./(1 + exp(-2*rx_cw./(N0/2)));
P0 = 1 - P1;

% Init
K0 = zeros(M, N);
K1 = zeros(M, N);
R0 = zeros(M, N);
R1 = zeros(M, N);
Q0 = H.*repmat(P0', M, 1);
Q1 = H.*repmat(P1', M, 1);
uhat = zeros(N,1);

for n = 1:Niter

    % loop over each check node
    for ii = 1:M

        % Find non-zeros in the row, i.e., find the connected variable nodes
        c1 = find(H(ii, :));

        for k = 1:length(c1)

            % Get column products of dR\c1(l) in step 2
            % Q0,Q1 corresponds to p_{j-i}(b)
            dR = 1.0;
            for l = 1:length(c1)
                if l~= k
                    dR = dR*(Q0(ii, c1(l)) - Q1(ii, c1(l)));
                end
            end % for l
            
            % the p_{i-j}(0) and p_{i-j}(1) in step 2
            R0(ii, c1(k)) = (1 + dR)/2;
            R1(ii, c1(k)) = (1 - dR)/2;

        end % for k

    end % for ii

    % loop over each variable nodes
    for jj = 1:N

        % Find non-zeros in the column, i.e., find the connected check nodes
        r1 = find(H(:, jj));

        for k = 1:length(r1)

            % Get row products of prodOfrij\ri(l)
            prodOfrij0 = 1;
            prodOfrij1 = 1;
            for l = 1:length(r1)
                if l~= k
                    prodOfrij0 = prodOfrij0*R0(r1(l), jj);
                    prodOfrij1 = prodOfrij1*R1(r1(l), jj);
                end
            end % for l

            % Update constants
            K0(r1(k), jj) = P0(jj)*prodOfrij0;
            K1(r1(k), jj) = P1(jj)*prodOfrij1;

            % Update Q0 and Q1, ccorresponding to p_{j-i}(b) in step 3
            Q0(r1(k), jj) = K0(r1(k), jj)./(K0(r1(k), jj) + K1(r1(k), jj));
            Q1(r1(k), jj) = K1(r1(k), jj)./(K0(r1(k), jj) + K1(r1(k), jj));

        end % for k

        % Update constants, this is step 4
        Ki0 = P0(jj)*prod(R0(r1, jj));
        Ki1 = P1(jj)*prod(R1(r1, jj));

        % Get Qj
        Qi0 = Ki0/(Ki0 + Ki1);
        Qi1 = Ki1/(Ki0 + Ki1);

        % Decode Qj, step 5
        if Qi1 > Qi0
            uhat(jj) = 1;
        else
            uhat(jj) = 0;
        end

    end % for nx

    % check if the decoding is already successful
    symdrome = mod(H*uhat, 2);
    if isempty(find(symdrome, 1))
        break
    end
    
end % for n
