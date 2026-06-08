function [sub_nodes, sub_adj_list, sub_E, sub_edge_weights, target_local_idx] = ...
    khop_subgraph(target_node, k, adj_list, E, edge_weights)
% KHOP_SUBGRAPH  Extract k-hop neighborhood subgraph around a target node.
%
%   For a k-layer GINEConv GNN, node target_node's output depends only on
%   its k-hop neighborhood. This function extracts that neighborhood and
%   returns a self-contained subgraph with remapped indices.
%
%   The subgraph approach is exact (not an approximation): running GINEConv
%   reachability on the subgraph produces identical output for target_node
%   as running on the full graph.
%
%   Reference: CORA GNN verification (TMLR), Corollary 17 (Subgraph Selection).
%              SCIP-MPNN (ICML 2024), k_hop_subgraph utility.
%
% Inputs:
%   target_node  - scalar, 1-indexed node to center on
%   k            - number of hops (= number of GINEConv layers in GNN)
%   adj_list     - (m x 2) edge list [src, dst], 1-indexed
%   E            - (m x E_in) edge feature matrix
%   edge_weights - (m x 1) aggregation weights
%
% Outputs:
%   sub_nodes        - (n_sub x 1) sorted original node indices in subgraph
%   sub_adj_list     - (m_sub x 2) remapped edge list (1:n_sub indexing)
%   sub_E            - (m_sub x E_in) edge features for subgraph edges
%   sub_edge_weights - (m_sub x 1) edge weights for subgraph edges
%   target_local_idx - target node's 1-indexed position in sub_nodes
%
% Author: Anne Tumlin
% Date: 03/11/2026

    if isempty(adj_list)
        error('khop_subgraph: adj_list is empty');
    end

    N = max(max(adj_list));

    % Build undirected neighbor lookup from directed edge list.
    % Power grid adj_lists store both (i,j) and (j,i), but we check both
    % columns to handle any directed graph correctly.
    neighbors = cell(N, 1);
    for e = 1:size(adj_list, 1)
        s = adj_list(e, 1);
        d = adj_list(e, 2);
        neighbors{s}(end+1) = d;
        neighbors{d}(end+1) = s;
    end

    % BFS for k hops from target_node
    visited = false(N, 1);
    visited(target_node) = true;
    frontier = target_node;

    for hop = 1:k
        new_frontier = [];
        for ni = 1:length(frontier)
            nbrs = neighbors{frontier(ni)};
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

    % Build global → local index map
    node_map = zeros(N, 1);
    node_map(sub_nodes) = 1:length(sub_nodes);
    target_local_idx = node_map(target_node);

    % Filter edges: keep only those where BOTH endpoints are in the subgraph
    edge_mask = visited(adj_list(:, 1)) & visited(adj_list(:, 2));
    sub_adj_raw = adj_list(edge_mask, :);

    % Remap endpoints to subgraph-local indices
    sub_adj_list = [node_map(sub_adj_raw(:, 1)), node_map(sub_adj_raw(:, 2))];

    if isa(E, 'GraphStar')
        sub_E = E.extractSubgraph(find(edge_mask));
    else
        sub_E = E(edge_mask, :);
    end

    if isempty(edge_weights)
        sub_edge_weights = [];
    else
        sub_edge_weights = edge_weights(edge_mask);
    end
end
