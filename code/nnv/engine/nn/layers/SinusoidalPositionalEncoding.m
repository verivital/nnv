classdef SinusoidalPositionalEncoding < handle
    % SinusoidalPositionalEncoding - Fixed sinusoidal positional encoding
    %
    % Implements the positional encoding from "Attention Is All You Need":
    %   PE(pos, 2i) = sin(pos / 10000^(2i/d_model))
    %   PE(pos, 2i+1) = cos(pos / 10000^(2i/d_model))
    %
    % This encoding is deterministic and fixed (not learned).
    % For classification tasks with fixed sequence length, this becomes
    % a constant addition to the input embeddings.
    %
    % Reference: https://arxiv.org/abs/1706.03762
    %
    % Author: NNV Team
    % Date: November 2025

    properties
        Name = 'sinusoidal_pe';
        NumInputs = 1;
        InputNames = {'in'};
        NumOutputs = 1;
        OutputNames = {'out'};

        % Parameters
        MaxSeqLen = 512;    % Maximum sequence length
        EmbedDim = 64;      % Embedding dimension (d_model)
        BaseFreq = 10000;   % Base frequency for sinusoidal encoding

        % Precomputed encoding matrix [MaxSeqLen x EmbedDim]
        PE = [];
    end

    methods

        function obj = SinusoidalPositionalEncoding(varargin)
            % Constructor
            % Usage:
            %   SinusoidalPositionalEncoding(embed_dim)
            %   SinusoidalPositionalEncoding(embed_dim, max_seq_len)
            %   SinusoidalPositionalEncoding(embed_dim, max_seq_len, base_freq)
            %   SinusoidalPositionalEncoding(embed_dim, max_seq_len, base_freq, name)

            switch nargin
                case 0
                    % Default
                case 1
                    obj.EmbedDim = varargin{1};
                case 2
                    obj.EmbedDim = varargin{1};
                    obj.MaxSeqLen = varargin{2};
                case 3
                    obj.EmbedDim = varargin{1};
                    obj.MaxSeqLen = varargin{2};
                    obj.BaseFreq = varargin{3};
                case 4
                    obj.EmbedDim = varargin{1};
                    obj.MaxSeqLen = varargin{2};
                    obj.BaseFreq = varargin{3};
                    obj.Name = varargin{4};
                otherwise
                    error('Invalid number of arguments');
            end

            % Precompute the positional encoding matrix
            obj.PE = obj.compute_encoding();
        end

        function PE = compute_encoding(obj)
            % Compute the sinusoidal positional encoding matrix
            % Returns: PE [MaxSeqLen x EmbedDim]

            PE = zeros(obj.MaxSeqLen, obj.EmbedDim);

            position = (0:obj.MaxSeqLen-1)';  % [MaxSeqLen x 1]

            % Compute division terms: 10000^(2i/d_model) for i = 0, 1, ...
            div_term = exp((0:2:obj.EmbedDim-1) * (-log(obj.BaseFreq) / obj.EmbedDim));

            % Apply sin to even indices, cos to odd indices
            PE(:, 1:2:end) = sin(position * div_term);
            PE(:, 2:2:end) = cos(position * div_term(1:floor(obj.EmbedDim/2)));

            % Handle odd embedding dimension
            if mod(obj.EmbedDim, 2) == 1
                % Last column is sin only
                PE(:, end) = sin(position * div_term(end));
            end
        end

    end

    methods  % Evaluation

        function y = evaluate(obj, x, positions)
            % Add positional encoding to input embeddings
            % @x: input embeddings [EmbedDim x SeqLen] or [EmbedDim x 1]
            % @positions: (optional) position indices, default is 0:SeqLen-1
            % @y: output with positional encoding added

            if nargin < 3
                % Default: positions are 0, 1, 2, ...
                if isvector(x)
                    seq_len = 1;
                    positions = 0;
                else
                    seq_len = size(x, 2);
                    positions = 0:seq_len-1;
                end
            end

            positions = positions + 1;  % MATLAB 1-indexed

            % Get positional encodings for these positions
            pe_subset = obj.PE(positions, :)';  % [EmbedDim x SeqLen]

            % Add to input
            if isvector(x)
                y = x(:) + pe_subset(:);
            else
                y = x + pe_subset;
            end
        end

        function pe = get_encoding(obj, positions)
            % Get positional encoding for specific positions
            % @positions: position indices (0-indexed)
            % @pe: positional encoding vectors [EmbedDim x length(positions)]

            positions = positions + 1;  % MATLAB 1-indexed
            pe = obj.PE(positions, :)';
        end

    end

    methods  % Reachability

        function R = reach(obj, varargin)
            % Main reachability method
            % Since PE is a fixed constant addition, reachability is exact

            switch nargin
                case 2
                    I = varargin{1};
                    method = 'approx-star';
                    positions = [];
                case 3
                    I = varargin{1};
                    method = varargin{2};
                    positions = [];
                case 4
                    I = varargin{1};
                    method = varargin{2};
                    positions = varargin{3};
                otherwise
                    I = varargin{1};
                    method = varargin{2};
                    positions = [];
            end

            if isa(I, 'ImageStar')
                R = obj.reach_imagestar(I, positions);
            elseif isa(I, 'Star')
                R = obj.reach_star(I, positions);
            else
                error('Input must be Star or ImageStar');
            end
        end

        function R = reach_star(obj, I, positions)
            % Reachability for Star input
            % Adding a constant vector to a Star is exact (just shift the center)

            if ~isa(I, 'Star')
                error('Input must be a Star');
            end

            n = I.dim;

            % Determine positions
            if isempty(positions)
                % Assume single position (position 0)
                positions = 0;
            end

            % Get positional encoding
            pe = obj.get_encoding(positions);
            pe = pe(:);  % Flatten to match Star dimension

            % Ensure dimensions match
            if length(pe) ~= n
                % If dimensions don't match, we might be dealing with
                % a different input format. Try to broadcast or error.
                if n == obj.EmbedDim
                    pe = pe(1:n);
                else
                    error('Dimension mismatch: Star dim=%d, PE dim=%d', n, length(pe));
                end
            end

            % Adding constant to Star: just shift the center
            % S = {c + V*alpha : C*alpha <= d}
            % S + pe = {(c + pe) + V*alpha : C*alpha <= d}

            new_V = I.V;
            new_V(:, 1) = new_V(:, 1) + pe;  % Add to center

            R = Star(new_V, I.C, I.d, I.predicate_lb, I.predicate_ub);
        end

        function R = reach_imagestar(obj, I, positions)
            % Reachability for ImageStar input

            if ~isa(I, 'ImageStar')
                error('Input must be an ImageStar');
            end

            % Convert to Star, process, convert back
            S = I.toStar;
            R_star = obj.reach_star(S, positions);

            % Convert back to ImageStar
            R = R_star.toImageStar(I.height, I.width, I.numChannel);
        end

        function R = reach_star_single_input(obj, I, method, ~, ~, ~)
            % Compatibility method
            R = obj.reach_star(I, []);
        end

    end

    methods(Static)

        function L = parse(layer)
            % Parse - not typically used as this is a custom layer
            error('SinusoidalPositionalEncoding parsing not implemented');
        end

    end

    methods  % Helper methods

        function obj = toGPU(obj)
            obj.PE = gpuArray(obj.PE);
        end

        function obj = changeParamsPrecision(obj, precision)
            if strcmp(precision, 'single')
                obj.PE = single(obj.PE);
            elseif strcmp(precision, 'double')
                obj.PE = double(obj.PE);
            end
        end

    end

end
