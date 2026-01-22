classdef LearnedPositionalEmbedding < handle
    % LearnedPositionalEmbedding - Learned positional embedding layer
    %
    % Learned position embeddings map each position index to a learnable
    % embedding vector. Unlike sinusoidal encodings, these are trained
    % with the model.
    %
    % Usage:
    %   layer = LearnedPositionalEmbedding(embed_dim, max_seq_len)
    %   layer = LearnedPositionalEmbedding(PE_matrix)
    %   layer = LearnedPositionalEmbedding(PE_matrix, name)
    %
    % Properties:
    %   PE - Positional embedding matrix [MaxSeqLen x EmbedDim]
    %   EmbedDim - Embedding dimension
    %   MaxSeqLen - Maximum sequence length
    %
    % Author: NNV Team
    % Date: November 2025

    properties
        Name = 'learned_pe';
        PE              % Positional embedding matrix [MaxSeqLen x EmbedDim]
        EmbedDim        % Embedding dimension
        MaxSeqLen       % Maximum sequence length
        NumInputs = 1;
        NumOutputs = 1;
        InputNames = {'in'};
        OutputNames = {'out'};
    end

    methods

        function obj = LearnedPositionalEmbedding(varargin)
            % Constructor
            %
            % Forms:
            %   LearnedPositionalEmbedding(embed_dim, max_seq_len)
            %   LearnedPositionalEmbedding(PE_matrix)
            %   LearnedPositionalEmbedding(PE_matrix, name)
            %   LearnedPositionalEmbedding(embed_dim, max_seq_len, name)

            if nargin == 0
                return;
            end

            if nargin == 1
                % Single argument: PE matrix
                obj.PE = varargin{1};
                obj.MaxSeqLen = size(obj.PE, 1);
                obj.EmbedDim = size(obj.PE, 2);

            elseif nargin == 2
                if ischar(varargin{2}) || isstring(varargin{2})
                    % (PE_matrix, name)
                    obj.PE = varargin{1};
                    obj.MaxSeqLen = size(obj.PE, 1);
                    obj.EmbedDim = size(obj.PE, 2);
                    obj.Name = char(varargin{2});
                else
                    % (embed_dim, max_seq_len) - initialize randomly
                    obj.EmbedDim = varargin{1};
                    obj.MaxSeqLen = varargin{2};
                    % Initialize with small random values (Xavier-like)
                    obj.PE = randn(obj.MaxSeqLen, obj.EmbedDim) * sqrt(2 / (obj.MaxSeqLen + obj.EmbedDim));
                end

            elseif nargin == 3
                % (embed_dim, max_seq_len, name)
                obj.EmbedDim = varargin{1};
                obj.MaxSeqLen = varargin{2};
                obj.Name = char(varargin{3});
                obj.PE = randn(obj.MaxSeqLen, obj.EmbedDim) * sqrt(2 / (obj.MaxSeqLen + obj.EmbedDim));
            end
        end

        function pe = get_encoding(obj, position)
            % Get positional encoding for a single position
            % @position: 0-indexed position
            % @pe: embedding vector [EmbedDim x 1]

            if position < 0 || position >= obj.MaxSeqLen
                error('Position %d out of range [0, %d)', position, obj.MaxSeqLen);
            end

            % Convert to 1-indexed
            pe = obj.PE(position + 1, :)';
        end

        function y = evaluate(obj, x, positions)
            % Add positional embeddings to input
            % @x: input embeddings [EmbedDim x SeqLen]
            % @positions: 0-indexed position indices [1 x SeqLen]
            % @y: output with PE added [EmbedDim x SeqLen]

            if nargin < 3
                % Default: positions 0, 1, 2, ...
                seq_len = size(x, 2);
                positions = 0:(seq_len - 1);
            end

            if any(positions < 0) || any(positions >= obj.MaxSeqLen)
                error('Positions out of range [0, %d)', obj.MaxSeqLen);
            end

            % Get positional embeddings: [EmbedDim x SeqLen]
            pe = obj.PE(positions + 1, :)';  % Convert to 1-indexed

            % Add to input
            y = x + pe;
        end

        function S_out = reach(obj, S_in, method, position)
            % Reachability analysis
            % Adding a constant (positional encoding) to a Star just shifts the center
            %
            % @S_in: input Star set
            % @method: reachability method (ignored, exact for constant addition)
            % @position: position index (0-indexed) or array of positions
            % @S_out: output Star set

            if nargin < 3
                method = 'approx-star';
            end
            if nargin < 4
                % Default: single position 0
                position = 0;
            end

            if isa(S_in, 'Star')
                S_out = obj.reach_star(S_in, position);
            elseif isa(S_in, 'ImageStar')
                S_out = obj.reach_imagestar(S_in, position);
            else
                error('Input must be Star or ImageStar');
            end
        end

        function S_out = reach_star(obj, S_in, position)
            % Reachability for Star set
            % Adding constant to Star: shift the center

            if isscalar(position)
                % Single position
                pe = obj.get_encoding(position);

                % Shift center by positional encoding
                new_c = S_in.V(:, 1) + pe;
                new_V = [new_c, S_in.V(:, 2:end)];

                S_out = Star(new_V, S_in.C, S_in.d, S_in.predicate_lb, S_in.predicate_ub);
            else
                % Multiple positions - flatten and concatenate
                seq_len = length(position);
                pe_all = obj.PE(position + 1, :)';  % [EmbedDim x SeqLen]
                pe_flat = pe_all(:);  % Flatten

                % Shift center
                new_c = S_in.V(:, 1) + pe_flat;
                new_V = [new_c, S_in.V(:, 2:end)];

                S_out = Star(new_V, S_in.C, S_in.d, S_in.predicate_lb, S_in.predicate_ub);
            end
        end

        function IS_out = reach_imagestar(obj, IS_in, position)
            % Reachability for ImageStar

            pe = obj.get_encoding(position);

            % Reshape PE to match ImageStar dimensions
            pe_img = reshape(pe, IS_in.height, IS_in.width, IS_in.numChannel);

            % Shift center
            new_V = IS_in.V;
            new_V(:, :, :, 1) = IS_in.V(:, :, :, 1) + pe_img;

            IS_out = ImageStar(new_V, IS_in.C, IS_in.d, IS_in.pred_lb, IS_in.pred_ub);
        end

        function obj = set_embeddings(obj, new_PE)
            % Set new positional embeddings
            % @new_PE: new embedding matrix [MaxSeqLen x EmbedDim]

            obj.PE = new_PE;
            obj.MaxSeqLen = size(new_PE, 1);
            obj.EmbedDim = size(new_PE, 2);
        end

        function obj = initialize_random(obj, embed_dim, max_seq_len)
            % Initialize with random embeddings

            obj.EmbedDim = embed_dim;
            obj.MaxSeqLen = max_seq_len;
            obj.PE = randn(max_seq_len, embed_dim) * sqrt(2 / (max_seq_len + embed_dim));
        end

        function obj = initialize_zeros(obj, embed_dim, max_seq_len)
            % Initialize with zeros (common for learned PE)

            obj.EmbedDim = embed_dim;
            obj.MaxSeqLen = max_seq_len;
            obj.PE = zeros(max_seq_len, embed_dim);
        end

        function sim = position_similarity(obj, pos1, pos2)
            % Compute cosine similarity between two position embeddings

            pe1 = obj.get_encoding(pos1);
            pe2 = obj.get_encoding(pos2);

            sim = dot(pe1, pe2) / (norm(pe1) * norm(pe2));
        end

        function info = get_info(obj)
            % Get layer information

            info = struct();
            info.Name = obj.Name;
            info.Type = 'LearnedPositionalEmbedding';
            info.EmbedDim = obj.EmbedDim;
            info.MaxSeqLen = obj.MaxSeqLen;
            info.NumParameters = obj.MaxSeqLen * obj.EmbedDim;
        end

    end

    methods (Static)

        function layer = from_pretrained(weights, name)
            % Create layer from pretrained weights
            % @weights: PE matrix [MaxSeqLen x EmbedDim]
            % @name: optional layer name

            if nargin < 2
                name = 'learned_pe';
            end

            layer = LearnedPositionalEmbedding(weights, name);
        end

    end

end
