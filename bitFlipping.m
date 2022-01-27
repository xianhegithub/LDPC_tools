function uhat = bitFlipping(rx_cw, H, Niter)
% Hard-decision/bit flipping sum product algorithm LDPC decoder
%
% Inputs
%  rx_cw     : Received signal vector (column vector)
%  H         : LDPC matrix
%  Niter     : num of iterations
%
% Outputs
%  uhat      : Decoded vector (0/1)

[M, N] = size(H);

% Prior hard-decision
ci = 0.5*(sign(rx_cw') + 1);

% Init
% R is the check-node matrix
R = zeros(M, N);
uhat = zeros(N,1);

% Asscociate the ci matrix with non-zero elements of H
Q = H.*repmat(ci, M, 1);

for n = 1:Niter


    % index over rows
    for ii = 1:M

        % Find non-zeros in the rwo ii
        c1 = find(H(ii, :));

        % Get the summation of Q\c1(k)
        for k = 1:length(c1)

            R(ii, c1(k)) = mod(sum(Q(ii, c1)) + Q(ii, c1(k)), 2);

        end % for k

    end % for i

    % index over columns
    for jj = 1:N

        % Find non-zero in the column jj
        r1 = find(H(:, jj));

        % Number of 1s in a column
        numOfOnes = length(find(R(r1, jj)));

        for k = 1:length(r1)

            % Update Q, set '1' for majority of 1s else '0', excluding r1(k)
            if numOfOnes + ci(jj) >= length(r1) - numOfOnes + R(r1(k), jj)
                Q(r1(k), jj) = 1;
            else
                Q(r1(k), jj) = 0;
            end

        end % for k

        % Bit decoding
        if numOfOnes + ci(jj) >= length(r1) - numOfOnes
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