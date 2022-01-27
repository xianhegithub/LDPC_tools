function H = PEG_ACE(N, M, Dv, DoYouWantACE)
% This function construct a LDPC parity-check matrix using the progressive
% edge growth algorithm
% Random PEG. If you want add ACE improvement, just set 'DoYouWantACE = 1'.
% Ref. [1] Regular and Irregular Progressive Edge-Growth Tanner Graphs --
% Xiao-Yu Hu et al.
% data structure from [2] https://github.com/YuYongRun/LDPC
% input: Dv degrees of each variable node, which can be randomly generated according to
% the varaiable degree distribution

H = zeros(M, N, 'uint8');                                                  % Unsigned integer 8 is enough.
% The following four arrays are in fact not necessary. However, if you use it, the implementation speed will significantly improve. It is some sort
% of sparse representation of H.
vn_cxn = zeros(max(Dv), N);
vn_deg_curr = ones(N, 1);
cn_cxn = zeros(M, 2 * ceil(sum(Dv)/M));
cn_deg_curr = ones(M, 1);

tree = zeros(N + M, 1, 'uint16');                                          % N + M is the upper bound of the number of elements in the tree.
parent = zeros(N + M, 1, 'uint16');                                        % Record the relationships in the tree.
CN_location = zeros(M, 1);                                                 % Record the relationships in the tree, this is only useful for ACE
VN_location = zeros(N, 1);                                                 % Record the relationships in the tree, this is only useful for ACE
%'tree' and 'parent' share the same counter.
marker = zeros(2 * M - 2, 2);                                              % 2M - 2 is the upper bound of the number of elements in the marker.
for curr_vn = 1 : N
    if mod(curr_vn, ceil(N/20)) == 1                                       % this is just a progress bar.
        disp(['PEG LDPC under Construction. Now ' num2str(curr_vn/N*100) '%.']);
    end
    for k = 1 : Dv(curr_vn)%for each k, we add an '1' in H.
        if k == 1                                                          % get the 0-th level for k = 1
            
            row_weight = sum(H, 2);
            cand_cn = find(row_weight == min(row_weight));
            idx = unidrnd(length(cand_cn));                                % random selection.
            sel_cn_idx = cand_cn(idx);
            H(sel_cn_idx, curr_vn) = 1;                                    % add an edge.
            vn_cxn(vn_deg_curr(curr_vn), curr_vn) = sel_cn_idx;            % record the H_col_rec_idx(curr_vn) th connection of the current variable node is the selected check node
            vn_deg_curr(curr_vn) = vn_deg_curr(curr_vn) + 1;               % increase the current degree of the current variable node by 1
            cn_cxn(sel_cn_idx, cn_deg_curr(sel_cn_idx)) = curr_vn;         % record that the H_row_rec_idx(sel_cn_idx) th connection of the selected check node is the current variable node
            cn_deg_curr(sel_cn_idx) = cn_deg_curr(sel_cn_idx) + 1;         % increase the current degree of the selected check node by 1
            
        else                                                               % get the 0-th level for k > 1
            level = 0;                                                     % refresh
            tree(1) = curr_vn;                                             % for every variable node, we build a new tree
            parent(1) = NaN;                                               % The Root Node has no parent.
            tree_idx = 1;                                                  % refresh
            VN_location(curr_vn) = tree_idx;
            tree_idx = tree_idx + 1;
            marker(1, :) = [1 1];                                          % on the first level, variable node starts at idx 1, ends at idx 1
            cn_level = find(H(:, curr_vn));                                % return the index of check nodes currently connected to curr_vn
            marker(2, 1) = 2;                                              % the first check node starts at 2 in tree vector
            for i = 1 : length(cn_level)
                tree(tree_idx) = cn_level(i);
                parent(tree_idx) = curr_vn;
                CN_location(cn_level(i)) = tree_idx;
                tree_idx = tree_idx + 1;
            end
            marker(2, 2) = marker(2, 1) + length(cn_level) - 1;            % where the check nodes in this level ends in tree vector
            marker(3, 1) = marker(2, 2) + 1;
            marker(3, 2) = marker(2, 2);
        
            level = level + 1;
            while(1)                                                       % this is where you expand the subgraph from variable node s_j
                prev_cn = zeros(M, 1);                                     % in every iteration, prev_cn/vn is re-initiated?
                prev_vn = zeros(N, 1);
                for k1 = 0 : level - 1
                    for k2 = marker(2 * k1 + 1, 1) : marker(2 * k1 + 1, 2) % loop through all variable nodes on previous levels and mark it '1' in the vector prev_vn
                        prev_vn(tree(k2)) = 1;
                    end
                    for k2 = marker(2 * k1 + 2, 1) : marker(2 * k1 + 2, 2) % loop through all check nodes on previous levels and mark it '1' in the vector prev_cn
                        prev_cn(tree(k2)) = 1;
                    end
                end
                %Add VNs in this level to the tree
                cannot_grow_flag = 1;
                for i = marker(2 * level, 1) : marker(2 * level, 2)        % loop through the check nodes on the level
                    for j = 1 : cn_deg_curr(tree(i)) - 1                   % for each check node, loop through all its current connections and add the variable nodes
                        target_vn = cn_cxn(tree(i), j);
                        if ~prev_vn(target_vn)                             % add the target_vn to the tree, to the prev_vn list.
                            prev_vn(target_vn) = 1;
                            tree(tree_idx) = target_vn;
                            parent(tree_idx) = tree(i);
                            VN_location(target_vn) = tree_idx;
                            tree_idx = tree_idx + 1;
                            marker(2 * level + 1, 2) = marker(2 * level + 1, 2) + 1;
                            cannot_grow_flag = 0;                          % if a variable node is added to the tree, in most cases, it means we should expand the subgraph further along it
                        end
                    end
                end 
                if cannot_grow_flag && sum(prev_cn) < M                    % This level is empty, i.e., no VNs is connected to the check nodes on the previous level
                    prev_cn_cmpl = find(mod(prev_cn + 1, 2));              % find the zeroes in prev_cn, prev_cn_cmpl is the list of check nodes that are not in the tree.
                    row_weight = sum(H(prev_cn_cmpl, :), 2);               % find the check nodes with minimum degree in prev_cn_cmpl, if multiple, then random select one
                    cand_cn = find(row_weight == min(row_weight));
                    idx = unidrnd(length(cand_cn));                        % Random selection.
                    sel_cn_idx = cand_cn(idx);
                    H(prev_cn_cmpl(sel_cn_idx), curr_vn) = 1;              % pick the check node from the complementary of N^l_{s_j}
                    vn_cxn(vn_deg_curr(curr_vn), curr_vn) = prev_cn_cmpl(sel_cn_idx); 
                    vn_deg_curr(curr_vn) = vn_deg_curr(curr_vn) + 1;
                    cn_cxn(prev_cn_cmpl(sel_cn_idx), cn_deg_curr(prev_cn_cmpl(sel_cn_idx))) = curr_vn;
                    cn_deg_curr(prev_cn_cmpl(sel_cn_idx)) = cn_deg_curr(prev_cn_cmpl(sel_cn_idx)) + 1;
                    break;                                                 % break out of while(1). This is the first condition where subgraph stops growing
                end
                marker(2 * level + 2, 1) = marker(2 * level + 1, 2) + 1;
                marker(2 * level + 2, 2) = marker(2 * level + 1, 2);
                %Add CNs
                now_cn = prev_cn;
                for i = marker(2 * level + 1, 1) : marker(2 * level + 1, 2)% loop through the variable nodes on the level
                    for j = 1 : vn_deg_curr(tree(i)) - 1                   % for each variable node, loop through all its current connections and add the check nodes
                        target_cn = vn_cxn(j, tree(i));
                        if ~now_cn(target_cn)
                            now_cn(target_cn) = 1;
                            tree(tree_idx) = target_cn;
                            parent(tree_idx) = tree(i);
                            CN_location(target_cn) = tree_idx;
                            tree_idx = tree_idx + 1;
                            marker(2 * level + 2, 2) = marker(2 * level + 2, 2) + 1;
                        end
                    end
                end
% Here, it seems that we should decide whether there is no new CN being added into the tree. If so, the tree stops growing, and
% we should add an '1' into the H matrix. However, this operation is not necessary. The deliberately designed array 'marker' will
% automatically guide the PEG process to the case where there is no new VN being added in to the tree in the next level, i.e.,
% we only need to see whether there is no new VNs being added.
                if sum(now_cn) == M                                        % ALL CNs are included. Then we turn back to the lower level, i.e., L - 1, to get the target CN.
                    prev_cn_cmpl = find(mod(prev_cn + 1, 2));
                    row_weight = sum(H(prev_cn_cmpl, :), 2);
                    cand_cn = find(row_weight == min(row_weight));
                    sel_cn_idx = NaN;
                    if DoYouWantACE
                        max_ACE = -realmax;
                        for i_ace = 1 : length(cand_cn)
                            involved_vn = zeros(level + 1, 1);
                            involved_cn = zeros(level + 1, 1);
                            involved_cn(1) = prev_cn_cmpl(cand_cn(i_ace));
                            involved_vn(1) = parent(CN_location(involved_cn(1)));
                            for i_back = 2 : 2 * level + 1
                                if mod(i_back, 2) == 0
                                    involved_cn(i_back/2 + 1) = parent(VN_location(involved_vn(i_back/2)));
                                else
                                    involved_vn((i_back - 1)/2 + 1) = parent(CN_location(involved_cn((i_back - 1)/2 + 1)));
                                end
                            end
                            ACE = sum(sum(H(:, involved_vn))) - 2 * (level + 1);
                            if ACE > max_ACE
                                max_ACE = ACE;
                                sel_cn_idx = prev_cn_cmpl(cand_cn(i_ace));
                            end
                        end
                    else
                        idx = unidrnd(length(cand_cn));                    % Random selection.
                        sel_cn_idx = prev_cn_cmpl(cand_cn(idx));
                    end
                    H(sel_cn_idx, curr_vn) = 1;
                    vn_cxn(vn_deg_curr(curr_vn), curr_vn) = sel_cn_idx;
                    vn_deg_curr(curr_vn) = vn_deg_curr(curr_vn) + 1;
                    cn_cxn(sel_cn_idx, cn_deg_curr(sel_cn_idx)) = curr_vn;
                    cn_deg_curr(sel_cn_idx) = cn_deg_curr(sel_cn_idx) + 1;
                    break;                                                 % break out of while(1). This is the second condition where subgraph stops growing
                end
                marker(2 * level + 3, 1) = marker(2 * level + 2, 2) + 1;
                marker(2 * level + 3, 2) = marker(2 * level + 2, 2);
                level = level + 1;
            end % end while(1)
        end  % end if k == 1
    end  % end for k = 1:Dv(curr_vn)
end  % end for curr_vn = 1:N
