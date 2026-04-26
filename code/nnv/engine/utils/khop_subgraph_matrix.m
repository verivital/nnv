function [sub_nodes, sub_A_norm, target_local_idx] = khop_subgraph_matrix(target_node, k, A_norm)
% KHOP_SUBGRAPH_MATRIX  Extract k-hop neighborhood subgraph for matrix-based GNNs.
%
%   For a k-layer GCN or SAGEConv GNN, node target_node's output depends
%   only on its k-hop neighborhood. This function extracts that neighborhood
%   and returns the corresponding submatrix of A_norm.
%
%   The subgraph approach is exact (not an approximation): the entries in
%   A_norm were computed with full-graph degrees, and the submatrix
%   preserves those normalization factors. Boundary nodes (at distance k)
%   may compute incorrect intermediate representations, but those are
%   never used in the target node's computation tree.
%
%   Reference: CORA GNN verification (TMLR), Corollary 17.
%
% Inputs:
%   target_node - scalar, 1-indexed node to center on
%   k           - number of hops (= number of message-passing layers)
%   A_norm      - (N x N) normalized adjacency matrix (GCN: D^{-1/2}(A+I)D^{-1/2},
%                 SAGEConv: binary adjacency with self-loops)
%
% Outputs:
%   sub_nodes        - (n_sub x 1) sorted original node indices in subgraph
%   sub_A_norm       - (n_sub x n_sub) submatrix of A_norm for subgraph nodes
%   target_local_idx - target node's 1-indexed position in sub_nodes
%
% Author: Anne Tumlin
% Date: 03/13/2026

    N = size(A_norm, 1);

    if target_node < 1 || target_node > N
        error('khop_subgraph_matrix: target_node %d out of range [1, %d]', target_node, N);
    end

    % Build neighbor list from adjacency sparsity pattern
    % Use A_norm ~= 0 to find connected pairs (works for both GCN and SAGEConv)
    if issparse(A_norm)
        A_pattern = A_norm;
    else
        A_pattern = sparse(A_norm);
    end

    % BFS for k hops from target_node
    visited = false(N, 1);
    visited(target_node) = true;
    frontier = target_node;

    for hop = 1:k
        new_frontier = [];
        for ni = 1:length(frontier)
            node = frontier(ni);
            % Find neighbors: nonzero entries in row and column
            nbrs_row = find(A_pattern(node, :));
            nbrs_col = find(A_pattern(:, node))';
            nbrs = unique([nbrs_row, nbrs_col]);
            for nb = nbrs
                if ~visited(nb)
                    visited(nb) = true;
                    new_frontier(end+1) = nb; %#ok<AGROW>
                end
            end
        end
        frontier = new_frontier;
        if isempty(frontier)
            break;
        end
    end

    sub_nodes = sort(find(visited));  % sorted original indices

    % Build global -> local index map
    node_map = zeros(N, 1);
    node_map(sub_nodes) = 1:length(sub_nodes);
    target_local_idx = node_map(target_node);

    % Extract submatrix (preserves full-graph normalization)
    sub_A_norm = full(A_norm(sub_nodes, sub_nodes));
end
