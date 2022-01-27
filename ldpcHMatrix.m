function H = ldpcHMatrix(num_r, num_c, method, no4cycle, checksPerCol)
% Create rate 1/2 LDPC adjacency matrix
%
% Inputs
%  num_r        : num of rows M
%  num_c        : num of column N
%  method   : technique for distributing non-zero element
%             0 evencol : For each column, place 1s uniformly at random
%             1 evenboth: For each column and row, place 1s uniformly at random
%  no4cyle  : if 1, then eliminate length-4 cycles, otherwise do nothing
%  checksPerCol: num of ones per column
%
% Output
%  H        : LDPC adjacency matrix

% num of ones per row (M/N ratio/code rate must be 0.5)

% References
% [1] http://www.cs.toronto.edu/~radford/ftp/LDPC-2006-02-08/pchk.html (main)
% [2] http://arun-10.tripod.com/ldpc/generate.html

if num_r/num_c ~= 0.5
    fprintf('Code rate must be 1/2 for now\n');
end
% only used for evenboth (when method ==1)
checksPerRow = ceil(num_c/num_r*checksPerCol);    % not really clear wether it should be ceil or floor.

% [1] step 1, Create a preliminary parity check matrix by one of the methods.
onesInCol = zeros(num_r,num_c);
% Distribute 1's uniformly at random within column
for nn = 1:num_c
    onesInCol(:, nn) = randperm(num_r);
end
% Create non zero elements (1s) index
r_idx = reshape(onesInCol(1:checksPerCol, :), num_c*checksPerCol, 1);
tmp = repmat((1:num_c), checksPerCol, 1);
c_idx = reshape(tmp, num_c*checksPerCol, 1);

if method == 0      % evencol method
    % Create sparse matrix H
    allOnes = ones(num_c*checksPerCol,1);
    H = full(sparse(r_idx, c_idx, allOnes, num_r, num_c));
    
elseif method == 1  % evenboth method
    % Make the number of 1s between rows as uniform as possible     
    % Order row index
    [~, sort_r_idx] = sort(r_idx);
    % Order column index based on row index
    c_idx_sort = zeros(num_c*checksPerCol,1);
    for nn = 1:num_c*checksPerCol
        c_idx_sort(nn) = c_idx(sort_r_idx(nn));
    end
    
    % Create new row index with uniform weight
    tmp = repmat(1:num_r, checksPerRow, 1);
    r_idx = reshape(tmp, num_c*checksPerCol, 1);
    
    % Create  H
    allOnes = ones(num_c*checksPerCol,1);
    H = full(sparse(r_idx, c_idx_sort, allOnes, num_r, num_c));
    
else
    fprintf('Unknown method\n');
    
end

% [1] step 2, Add 1s to the parity check matrix in order to avoid rows that
% have no 1s in them, or which have only one 1 in them.
% Check rows that have no 1 or only have one 1
% this is more likely to happen when using method 1 and num_c is large
for mm = 1:num_r

    n = randperm(num_c);
    % Add two 1s if row has no 1
    if isempty(find(r_idx == mm, 1)) 
        H(mm, n(1)) = 1;
        H(mm, n(2)) = 1;
    % Add one 1 if row has only one 1   
    elseif length(find(r_idx == mm)) == 1
        H(mm, n(1)) = 1;
    end
    
end % end mm = 1:num_r

% [1] step 3, if the preliminary parity check matrix constructed in step (1)
% had an even number of 1s in each column (checksPerCol is an even number),
% add further 1s to avoid the problem that this will cause the rows to add 
% to zero, and hence at least one check will be redundant.


% [1] step 4, eliminate situations where a pair of columns both have 1s in 
% a particular pair of rows, which correspond to cycles of length four in 
% the factor graph of the parity check matrix.
% if no4cycle, remove any length-4 cycle. This makes the code irregular.
if no4cycle == 1  

    for mm = 1:num_r

        % Look for pair of row - column
        for jx = (mm + 1):num_r
            w = H(mm, :)& H(jx, :);
            c1 = find(w);
            lc = length(c1);
            if lc > 1
                % If found, flip one 1 to 0 in the row with less number of 1s
                if length(find(H(mm, :))) < length(find(H(jx, :)))
                    % Repeat the process until only one column left 
                    for cc = 1:lc - 1
                        H(jx, c1(cc)) = 0;
                    end
                else
                    for cc = 1:lc - 1
                        H(mm, c1(cc)) = 0;
                    end
                end
            end
        end % for jx
    end % for mm
end % if no4cycle==1