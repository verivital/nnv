classdef RotaryPositionalEmbedding < handle
    % RotaryPositionalEmbedding - Rotary Position Embedding (RoPE) layer
    %
    % RoPE applies rotation matrices to pairs of dimensions based on position.
    % For each pair (2i, 2i+1), a 2D rotation is applied:
    %   x_rot[2i]   = x[2i] * cos(θ) - x[2i+1] * sin(θ)
    %   x_rot[2i+1] = x[2i] * sin(θ) + x[2i+1] * cos(θ)
    %
    % Where θ = position / base^(2i/dim)
    %
    % Key insight for verification: For fixed positions (classification),
    % RoPE is a LINEAR transformation (rotation matrix), making reachability
    % analysis straightforward.
    %
    % Usage:
    %   layer = RotaryPositionalEmbedding(dim)
    %   layer = RotaryPositionalEmbedding(dim, base)
    %   layer = RotaryPositionalEmbedding(dim, base, max_seq_len)
    %   y = layer.evaluate(x, positions)
    %   R = layer.reach(S_in, method, positions)
    %
    % Author: NNV Team
    % Date: November 2025
    %
    % Reference: RoFormer: Enhanced Transformer with Rotary Position Embedding
    %            https://arxiv.org/abs/2104.09864

    properties
        Name = 'rope';
        Dim             % Embedding dimension (must be even)
        Base = 10000;   % Base frequency (default 10000)
        MaxSeqLen = 512;% Maximum sequence length

        % Precomputed frequency components
        InvFreq         % Inverse frequencies [Dim/2 x 1]

        NumInputs = 1;
        NumOutputs = 1;
        InputNames = {'in'};
        OutputNames = {'out'};
    end

    methods

        function obj = RotaryPositionalEmbedding(varargin)
            % Constructor
            %
            % Forms:
            %   RotaryPositionalEmbedding(dim)
            %   RotaryPositionalEmbedding(dim, base)
            %   RotaryPositionalEmbedding(dim, base, max_seq_len)
            %   RotaryPositionalEmbedding(dim, base, max_seq_len, name)

            if nargin == 0
                return;
            end

            obj.Dim = varargin{1};

            if mod(obj.Dim, 2) ~= 0
                error('RoPE dimension must be even, got %d', obj.Dim);
            end

            if nargin >= 2
                obj.Base = varargin{2};
            end

            if nargin >= 3
                obj.MaxSeqLen = varargin{3};
            end

            if nargin >= 4
                obj.Name = char(varargin{4});
            end

            % Compute inverse frequencies
            obj.InvFreq = obj.compute_inv_freq();
        end

        function inv_freq = compute_inv_freq(obj)
            % Compute inverse frequency components
            % inv_freq[i] = 1 / base^(2i / dim)

            half_dim = obj.Dim / 2;
            inv_freq = 1.0 ./ (obj.Base .^ ((0:2:(obj.Dim-1))' / obj.Dim));
            % Result is [Dim/2 x 1]
        end

        function [cos_cached, sin_cached] = get_cos_sin(obj, positions)
            % Get cosine and sine values for given positions
            % @positions: array of position indices (0-indexed)
            % @cos_cached: [len(positions) x Dim/2]
            % @sin_cached: [len(positions) x Dim/2]

            positions = positions(:);  % Column vector

            % Compute angles: positions * inv_freq'
            % [num_pos x 1] * [1 x Dim/2] = [num_pos x Dim/2]
            freqs = positions * obj.InvFreq';

            cos_cached = cos(freqs);
            sin_cached = sin(freqs);
        end

        function y = evaluate(obj, x, positions)
            % Apply RoPE to input
            % @x: input tensor [Dim x SeqLen] or [Dim x 1]
            % @positions: 0-indexed position indices [1 x SeqLen]
            % @y: rotated output [Dim x SeqLen]

            if nargin < 3
                % Default positions: 0, 1, 2, ...
                seq_len = size(x, 2);
                positions = 0:(seq_len - 1);
            end

            seq_len = length(positions);
            half_dim = obj.Dim / 2;

            % Get cos and sin for positions
            [cos_pos, sin_pos] = obj.get_cos_sin(positions);
            % cos_pos, sin_pos: [SeqLen x Dim/2]

            y = zeros(size(x));

            for t = 1:seq_len
                % For each position, apply rotation to pairs
                cos_t = cos_pos(t, :)';  % [Dim/2 x 1]
                sin_t = sin_pos(t, :)';  % [Dim/2 x 1]

                % Split into even and odd indices
                x_even = x(1:2:end, t);  % [Dim/2 x 1]
                x_odd = x(2:2:end, t);   % [Dim/2 x 1]

                % Apply rotation
                y_even = x_even .* cos_t - x_odd .* sin_t;
                y_odd = x_even .* sin_t + x_odd .* cos_t;

                % Interleave back
                y(1:2:end, t) = y_even;
                y(2:2:end, t) = y_odd;
            end
        end

        function R = get_rotation_matrix(obj, position)
            % Get the rotation matrix for a single position
            % @position: 0-indexed position
            % @R: [Dim x Dim] block-diagonal rotation matrix

            [cos_pos, sin_pos] = obj.get_cos_sin(position);
            cos_pos = cos_pos(:);
            sin_pos = sin_pos(:);

            half_dim = obj.Dim / 2;
            R = zeros(obj.Dim, obj.Dim);

            for i = 1:half_dim
                idx_even = 2*i - 1;
                idx_odd = 2*i;

                % 2x2 rotation block
                R(idx_even, idx_even) = cos_pos(i);
                R(idx_even, idx_odd) = -sin_pos(i);
                R(idx_odd, idx_even) = sin_pos(i);
                R(idx_odd, idx_odd) = cos_pos(i);
            end
        end

        function S_out = reach(obj, S_in, method, positions)
            % Reachability analysis
            % For fixed positions, RoPE is a linear transformation!
            %
            % @S_in: input Star set
            % @method: reachability method (ignored - exact for linear)
            % @positions: position indices (scalar or array)
            % @S_out: output Star set

            if nargin < 3
                method = 'approx-star';
            end
            if nargin < 4
                positions = 0;  % Default: position 0
            end

            if isa(S_in, 'Star')
                S_out = obj.reach_star(S_in, positions);
            elseif isa(S_in, 'ImageStar')
                S_out = obj.reach_imagestar(S_in, positions);
            else
                error('Input must be Star or ImageStar');
            end
        end

        function S_out = reach_star(obj, S_in, positions)
            % Reachability for Star set
            % RoPE with fixed positions is LINEAR -> exact Star transformation

            if isscalar(positions)
                % Single position: apply rotation matrix
                R = obj.get_rotation_matrix(positions);

                % Linear transformation: S_out = R * S_in
                new_V = R * S_in.V;

                S_out = Star(new_V, S_in.C, S_in.d, S_in.predicate_lb, S_in.predicate_ub);
            else
                % Multiple positions: need to handle sequence
                % Input Star represents [Dim x SeqLen] flattened
                seq_len = length(positions);

                if S_in.dim ~= obj.Dim * seq_len
                    error('Star dimension (%d) must equal Dim * SeqLen (%d)', ...
                        S_in.dim, obj.Dim * seq_len);
                end

                % Build block-diagonal rotation matrix for all positions
                R_full = zeros(S_in.dim, S_in.dim);
                for t = 1:seq_len
                    R_t = obj.get_rotation_matrix(positions(t));
                    idx_start = (t-1) * obj.Dim + 1;
                    idx_end = t * obj.Dim;
                    R_full(idx_start:idx_end, idx_start:idx_end) = R_t;
                end

                % Apply transformation
                new_V = R_full * S_in.V;
                S_out = Star(new_V, S_in.C, S_in.d, S_in.predicate_lb, S_in.predicate_ub);
            end
        end

        function IS_out = reach_imagestar(obj, IS_in, positions)
            % Reachability for ImageStar (convert to Star, apply, convert back)

            % Flatten to Star
            S_in = IS_in.toStar();

            % Apply RoPE
            S_out = obj.reach_star(S_in, positions);

            % Convert back to ImageStar (if dimensions match)
            % This is a simplified version - may need adjustment for actual use
            IS_out = ImageStar(S_out.V, S_out.C, S_out.d, ...
                S_out.predicate_lb, S_out.predicate_ub);
        end

        function y = rotate_half(obj, x)
            % Rotate half operation (alternative RoPE implementation)
            % Swaps pairs and negates first half
            % Used in some implementations instead of interleaved rotation

            half = size(x, 1) / 2;
            y = [-x(half+1:end, :); x(1:half, :)];
        end

        function y = apply_rotary_emb_half(obj, x, positions)
            % Alternative RoPE implementation using rotate_half
            % Some models use this instead of interleaved

            [cos_pos, sin_pos] = obj.get_cos_sin(positions);

            % Expand cos/sin to full dimension by repeating
            cos_full = repmat(cos_pos', 2, 1);  % [Dim x SeqLen]
            sin_full = repmat(sin_pos', 2, 1);

            y = x .* cos_full + obj.rotate_half(x) .* sin_full;
        end

        function info = get_info(obj)
            % Get layer information

            info = struct();
            info.Name = obj.Name;
            info.Type = 'RotaryPositionalEmbedding';
            info.Dim = obj.Dim;
            info.Base = obj.Base;
            info.MaxSeqLen = obj.MaxSeqLen;
            info.NumParameters = 0;  % RoPE has no learnable parameters
        end

        function obj = set_dim(obj, dim)
            % Update dimension and recompute frequencies

            if mod(dim, 2) ~= 0
                error('RoPE dimension must be even');
            end
            obj.Dim = dim;
            obj.InvFreq = obj.compute_inv_freq();
        end

    end

    methods (Static)

        function test_rotation_properties()
            % Test that rotation matrices have expected properties

            rope = RotaryPositionalEmbedding(8);
            R = rope.get_rotation_matrix(5);

            % Rotation matrix should be orthogonal: R * R' = I
            err = norm(R * R' - eye(8), 'fro');
            assert(err < 1e-10, 'Rotation matrix not orthogonal');

            % Determinant should be 1
            assert(abs(det(R) - 1) < 1e-10, 'Determinant not 1');

            fprintf('Rotation matrix properties verified\n');
        end

    end

end
