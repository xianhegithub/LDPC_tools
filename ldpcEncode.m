function [checkBits, newH] = ldpcEncode(src, Hin, enc_strategy)
% Make parity check matrix from adjacency matrix H using LU
% decomposition, systematically encode the message blocks. 
%
% Inputs
%  src     : Binary source
%  H       : adjacency matrix
%  enc_strategy: Strategy for finding the next non-zero diagonal elements
%               0: First non-zero found by column search
%               1: Minimum number of non-zeros in later columns
%               2: Minimum product of:
%               - num of non-zeros its column minus 1
%               - num of non-zeros its row minus 1
% explanation of strategies can be found here
% http://www.cs.toronto.edu/~radford/ftp/LDPC-2006-02-08/sparse-LU.html
% Outputs
%  c       : Check bits
%  newH    : Rearrange H

% copy so as to not overwrite Hin
H = Hin; 
% Get the matric dimension
[num_r, num_c] = size(H);
% Set a new matrix F for LU decomposition
F = H;
% LU matrices
L = zeros(num_r, num_c - num_r);
U = zeros(num_r, num_c - num_r);

% Re-order the num_r x (num_c - num_r) submatrix
for mx = 1:num_r

    % strategy {0 = First; 1 = Mincol; 2 = Minprod}
    if enc_strategy == 0 % Create diagonal matrix using 'First' strategy
        % this strategy normally produces singular matrix, causing error

        % Find non-zero elements (1s) for the diagonal
        [r_idx, c_idx] = find(F(:, mx:end)~=0);

        % Find non-zero diagonal element candidates
        nz_r_idx = find(r_idx == c_idx);

        % Find the first non-zero column
        chosenCol = c_idx(nz_r_idx(1)) + (mx - 1);
        
    elseif enc_strategy== 1 % Create diag matrix using 'Mincol' strategy

        % Find non-zero elements (1s) for the diagonal
        [r_idx, c_idx] = find(F(:, mx:end)~=0);
        colWeight = sum(F(:, mx:end), 1);

        % Find non-zero diagonal element candidates
        [rowIndex, ~] = find(r_idx == mx);

        % Find the minimum column weight
        [~, ix] = min(colWeight(c_idx(rowIndex)));
        % Add offset to the chosen row index to match the dimension of the...
        % original matrix F
        chosenCol = c_idx(rowIndex(ix)) + (mx - 1);
        
    elseif enc_strategy == 2 % Create diag matrix using 'Minprod' strategy
        
        % Find non-zero elements of F
        [r_idx, c_idx] = find(F(:, mx:end)~=0);
        colWeight = sum(F(:, mx:end), 1) - 1;
        rowWeight = sum(F(mx, :), 2) - 1;

        % Find non-zero diagonal element candidates
        [rowIndex, ~] = find(r_idx == mx);

        % Find the minimum product
        [~, ix] = min(colWeight(c_idx(rowIndex))*rowWeight);
        % Add offset to the chosen row index to match the dimension of the
        % original matrix F
        chosenCol = c_idx(rowIndex(ix)) + (mx - 1);
        
    else
        
        fprintf("Unknown columns re-ordering strategy!\n");
        
    end % end if end_strategy 
    
    % Re-ordering columns of both H and F
    tmp1 = F(:, mx);
    tmp2 = H(:, mx);
    F(:, mx) = F(:, chosenCol);
    H(:, mx) = H(:, chosenCol);
    F(:, chosenCol) = tmp1;
    H(:, chosenCol) = tmp2;

    % Fill the LU matrices column by column
    L(mx:end, mx) = F(mx:end, mx);
    U(1:mx, mx) = F(1:mx, mx);
    
    % There will be no rows operation at the last row
    if mx < num_r

        % Find the later rows with non-zero elements in column mx
        [r2, ~] = find(F((mx + 1):end, mx)~=0);
        % Add current row to the later rows which have a 1 in column mx
        F((mx + r2), :) = mod(F((mx + r2), :) + repmat(F(mx, :), length(r2), 1), 2);

    end % if mx < num_r    
end  % for mx = 1:num_r

z = mod(H(:, (num_c - num_r) + 1:end)*src, 2);

% Parity check vector found by solving sparse LU
checkBits = mod(U\(L\z), 2);

% Return the rearrange H
newH = H;
