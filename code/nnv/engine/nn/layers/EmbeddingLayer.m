classdef EmbeddingLayer < handle
    % EmbeddingLayer - Word/Token embedding lookup layer
    %
    % Converts discrete token indices to dense embedding vectors.
    %   output = E(token_ids, :)
    %
    % For verification, perturbations can be applied in the embedding space
    % after the lookup (since token indices are discrete).
    %
    % Perturbation models supported:
    %   1. L-inf perturbation in embedding space
    %   2. Synonym substitution (discrete token perturbation)
    %
    % Author: NNV Team
    % Date: November 2025

    properties
        Name = 'embedding';
        NumInputs = 1;
        InputNames = {'in'};
        NumOutputs = 1;
        OutputNames = {'out'};

        % Embedding matrix [VocabSize x EmbedDim]
        E = [];

        % Dimensions
        VocabSize = 0;
        EmbedDim = 0;
    end

    methods

        function obj = EmbeddingLayer(varargin)
            % Constructor
            % Usage:
            %   EmbeddingLayer(embedding_matrix)
            %   EmbeddingLayer(vocab_size, embed_dim) - random init
            %   EmbeddingLayer(embedding_matrix, name)

            switch nargin
                case 0
                    % Default
                case 1
                    if isscalar(varargin{1})
                        error('Provide both vocab_size and embed_dim, or embedding matrix');
                    else
                        obj.E = varargin{1};
                        [obj.VocabSize, obj.EmbedDim] = size(obj.E);
                    end
                case 2
                    if isscalar(varargin{1}) && isscalar(varargin{2})
                        % Random initialization
                        obj.VocabSize = varargin{1};
                        obj.EmbedDim = varargin{2};
                        obj.E = randn(obj.VocabSize, obj.EmbedDim) * 0.02;
                    else
                        % Embedding matrix + name
                        obj.E = varargin{1};
                        obj.Name = varargin{2};
                        [obj.VocabSize, obj.EmbedDim] = size(obj.E);
                    end
                case 3
                    obj.VocabSize = varargin{1};
                    obj.EmbedDim = varargin{2};
                    obj.Name = varargin{3};
                    obj.E = randn(obj.VocabSize, obj.EmbedDim) * 0.02;
                otherwise
                    error('Invalid number of arguments');
            end
        end

    end

    methods  % Evaluation

        function y = evaluate(obj, token_ids)
            % Look up embeddings for token indices
            % @token_ids: token indices [SeqLen x 1] or [1 x SeqLen] (1-indexed)
            % @y: embeddings [EmbedDim x SeqLen]

            token_ids = token_ids(:)';  % Ensure row vector

            % Validate indices
            if any(token_ids < 1) || any(token_ids > obj.VocabSize)
                error('Token indices out of range [1, %d]', obj.VocabSize);
            end

            % Lookup
            y = obj.E(token_ids, :)';  % [EmbedDim x SeqLen]
        end

        function y = evaluate_single(obj, token_id)
            % Look up embedding for single token
            % @token_id: single token index (1-indexed)
            % @y: embedding vector [EmbedDim x 1]

            if token_id < 1 || token_id > obj.VocabSize
                error('Token index out of range [1, %d]', obj.VocabSize);
            end

            y = obj.E(token_id, :)';
        end

    end

    methods  % Reachability with perturbation models

        function R = reach(obj, varargin)
            % Main reachability method
            % Since token lookup is discrete, we need a perturbation model

            error('Use reach_with_epsilon or reach_with_synonyms for reachability');
        end

        function R = reach_with_epsilon(obj, token_ids, epsilon)
            % Reachability with L-inf epsilon-ball around embeddings
            % @token_ids: token indices (1-indexed)
            % @epsilon: perturbation bound in embedding space
            % @R: Star set representing all possible perturbed embeddings

            token_ids = token_ids(:)';
            seq_len = length(token_ids);

            % Get base embeddings
            base_emb = obj.evaluate(token_ids);  % [EmbedDim x SeqLen]
            base_emb = base_emb(:);  % Flatten

            % Create epsilon-ball around embeddings
            dim = length(base_emb);
            lb = base_emb - epsilon;
            ub = base_emb + epsilon;

            R = Star(lb, ub);
        end

        function R = reach_with_synonyms(obj, token_id, synonym_ids)
            % Reachability with discrete synonym substitution
            % @token_id: original token index
            % @synonym_ids: allowed synonym indices (including original)
            % @R: Union of embedding points (as Stars)

            if ~ismember(token_id, synonym_ids)
                synonym_ids = [token_id, synonym_ids];
            end

            % Get embeddings for all synonyms
            embeddings = obj.E(synonym_ids, :)';  % [EmbedDim x NumSynonyms]

            % Create convex hull as Star
            % The convex hull of discrete points
            n_syn = length(synonym_ids);

            if n_syn == 1
                % Single point
                emb = embeddings(:, 1);
                R = Star(emb, emb);
            else
                % Convex hull of synonym embeddings
                % Use interval bounds as overapproximation
                lb = min(embeddings, [], 2);
                ub = max(embeddings, [], 2);
                R = Star(lb, ub);
            end
        end

        function R = reach_star_single_input(obj, I, method, ~, ~, ~)
            % Compatibility method - not directly applicable
            error('EmbeddingLayer requires token indices. Use reach_with_epsilon.');
        end

    end

    methods(Static)

        function L = parse(layer)
            % Parse MATLAB/ONNX embedding layer

            if isprop(layer, 'Weights')
                E = layer.Weights;
            elseif isprop(layer, 'EmbeddingDimension') && isprop(layer, 'NumWords')
                % Create from dimensions
                E = randn(layer.NumWords, layer.EmbeddingDimension) * 0.02;
            else
                error('Cannot parse embedding layer');
            end

            if isprop(layer, 'Name')
                name = layer.Name;
            else
                name = 'embedding';
            end

            L = EmbeddingLayer(E, name);
        end

    end

    methods  % Helper methods

        function obj = set_embeddings(obj, E)
            % Set embedding matrix
            obj.E = E;
            [obj.VocabSize, obj.EmbedDim] = size(E);
        end

        function obj = toGPU(obj)
            obj.E = gpuArray(obj.E);
        end

        function obj = changeParamsPrecision(obj, precision)
            if strcmp(precision, 'single')
                obj.E = single(obj.E);
            elseif strcmp(precision, 'double')
                obj.E = double(obj.E);
            end
        end

        function emb = get_embedding(obj, token_id)
            % Get embedding for a single token
            emb = obj.E(token_id, :)';
        end

        function sim = cosine_similarity(obj, token_id1, token_id2)
            % Compute cosine similarity between two token embeddings
            e1 = obj.E(token_id1, :);
            e2 = obj.E(token_id2, :);
            sim = dot(e1, e2) / (norm(e1) * norm(e2));
        end

        function [ids, sims] = find_nearest(obj, token_id, k)
            % Find k nearest neighbors by cosine similarity
            % @token_id: query token
            % @k: number of neighbors
            % @ids: neighbor token indices
            % @sims: similarity scores

            query = obj.E(token_id, :);
            query_norm = norm(query);

            sims_all = zeros(obj.VocabSize, 1);
            for i = 1:obj.VocabSize
                if i ~= token_id
                    sims_all(i) = dot(query, obj.E(i, :)) / (query_norm * norm(obj.E(i, :)));
                else
                    sims_all(i) = -inf;  % Exclude self
                end
            end

            [sims_sorted, idx_sorted] = sort(sims_all, 'descend');
            ids = idx_sorted(1:k);
            sims = sims_sorted(1:k);
        end

    end

end
