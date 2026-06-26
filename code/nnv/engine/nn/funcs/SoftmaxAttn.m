classdef SoftmaxAttn
    % SoftmaxAttn - sound reachability primitives for softmax self-attention.
    %
    % This is the attention analogue of NNV's func classes (PosLin, LogSig, ...):
    % it holds the *sound* set-propagation math used by the transformer layers
    % (DynamicMatmulLayer, ScaledDotProductAttentionLayer) so the layer files stay
    % thin and the math is unit-testable in isolation.
    %
    % Ported from the n2v sound star-set ViT verifier (branch
    % feat/vit-vnncomp2023-sound), with every method validated in MATLAB by
    % Monte-Carlo containment against the linprog LP oracle
    % (soundness_test_utils.verify_star_containment). See
    %   tests/nn/attention/test_SoftmaxAttn_*.m
    % and the design notes in n2v docs/theory/sound-vit-reach.md.
    %
    % Conventions:
    %   - A "matrix Star" is a Star whose dim = prod(shape); its state vector is the
    %     column-major (MATLAB reshape) flattening of the matrix. All reshapes here
    %     use MATLAB's native column-major order, applied consistently to the center
    %     and every generator, so the predicate structure is preserved exactly.
    %   - "estimate" bounds = predicate-box over-approx (Star.estimateRanges, no LP,
    %     exact when there are no constraint rows). "lp" bounds = Star.getRanges (LP).
    %
    % Author: NNV Team (ViT track). Date: 2026.

    methods (Static)

        % ----- bounds helper -----------------------------------------------
        function [lo, hi] = starBounds(S, mode)
            % Sound per-dimension [lo,hi] of a Star, as column vectors.
            if nargin < 2, mode = 'estimate'; end
            if strcmpi(mode, 'lp')
                [lo, hi] = S.getRanges();          % LP-exact (honors C)
            else
                [lo, hi] = S.estimateRanges();     % predicate-box over-approx
            end
            lo = double(lo(:)); hi = double(hi(:));
        end

        % ----- Rump midpoint-radius interval matrix product ----------------
        function [cl, cu] = intervalMatMul(al, au, bl, bu)
            % Sound enclosure of { A*B : al<=A<=au, bl<=B<=bu }.
            % al,au are [m k]; bl,bu are [k n]; cl,cu are [m n].
            %
            % Rump's formula: with midpoints Am,Bm and radii Ar,Br,
            %   C = Am*Bm  +/-  (|Am|*Br + Ar*|Bm| + Ar*Br).
            % Each summand A(i,k)*B(k,j) deviates from Am*Bm by at most
            % |Am|*Br + Ar*|Bm| + Ar*Br (componentwise), and summing over the
            % contracted index k gives the radius matrix. Sound for operands of
            % arbitrary sign; no cross-correlation between the operands assumed.
            Am = (al + au) / 2;  Ar = (au - al) / 2;
            Bm = (bl + bu) / 2;  Br = (bu - bl) / 2;
            Cm = Am * Bm;
            Cr = abs(Am) * Br + Ar * abs(Bm) + Ar * Br;
            cl = Cm - Cr;
            cu = Cm + Cr;
        end

        % ----- set@set bilinear matmul (concretize / IBP-class) ------------
        function S = bilinearMatMulStar(X, Y, shX, shY, scale, mode)
            % Sound reach of S = scale * (Xmat * Ymat) where Xmat,Ymat are the
            % column-major reshapes of the matrix-Stars X (dim=prod(shX)) and
            % Y (dim=prod(shY)), shX=[m k], shY=[k n]. Returns a box Star of
            % dim m*n. Concretizes each operand to its sound per-dim range then
            % applies intervalMatMul -- sound, drops only operand cross-correlation.
            if nargin < 5 || isempty(scale), scale = 1.0; end
            if nargin < 6 || isempty(mode),  mode = 'estimate'; end
            if X.dim ~= prod(shX), error('SoftmaxAttn:shape','X.dim=%d ~= prod(shX)=%d', X.dim, prod(shX)); end
            if Y.dim ~= prod(shY), error('SoftmaxAttn:shape','Y.dim=%d ~= prod(shY)=%d', Y.dim, prod(shY)); end
            if shX(2) ~= shY(1), error('SoftmaxAttn:innerdim','inner dims differ %d vs %d', shX(2), shY(1)); end

            [xl, xu] = SoftmaxAttn.starBounds(X, mode);
            [yl, yu] = SoftmaxAttn.starBounds(Y, mode);
            Al = reshape(xl, shX); Au = reshape(xu, shX);
            Bl = reshape(yl, shY); Bu = reshape(yu, shY);
            [cl, cu] = SoftmaxAttn.intervalMatMul(Al, Au, Bl, Bu);
            % sign-safe scaling
            lo = min(scale .* cl, scale .* cu);
            hi = max(scale .* cl, scale .* cu);
            S = SoftmaxAttn.boxStar(reshape(lo, [], 1), reshape(hi, [], 1));
        end

        % ----- exact correlated row-softmax bound --------------------------
        function [a_lb, a_ub] = correlatedRowSoftmaxBounds(s_lo, s_hi)
            % Exact per-element range of A = softmax(S) over the KEY axis (dim 2)
            % given the elementwise logit box s_lo <= S <= s_hi (both [R x n]).
            %
            % softmax_j is increasing in its own logit (dA_j/dS_j = A_j(1-A_j) >= 0)
            % and decreasing in every other logit (dA_j/dS_k = -A_j A_k <= 0, k~=j),
            % so coordinate j attains its max at (s_lo off-j, s_hi at j) and its min
            % at (s_hi off-j, s_lo at j). These corners are evaluated exactly, giving
            % the TIGHTEST axis-aligned enclosure. Clamp to [0,1] to absorb FP.
            [R, n] = size(s_lo);
            a_lb = zeros(R, n);
            a_ub = zeros(R, n);
            for j = 1:n
                Su = s_lo;  Su(:, j) = s_hi(:, j);   % upper corner for column j
                a_ub(:, j) = SoftmaxAttn.softmaxColumn(Su, j);
                Sl = s_hi;  Sl(:, j) = s_lo(:, j);   % lower corner for column j
                a_lb(:, j) = SoftmaxAttn.softmaxColumn(Sl, j);
            end
            a_lb = min(max(a_lb, 0), 1);
            a_ub = min(max(a_ub, 0), 1);
        end

        function v = softmaxColumn(S, j)
            % numerically-stable row-softmax of [R x n] matrix S, return column j
            m = max(S, [], 2);
            E = exp(S - m);
            v = E(:, j) ./ sum(E, 2);
        end

        % ----- softmax-attention weights as a (box) Star -------------------
        function W = softmaxAttnStar(S_logits, shape, mode)
            % Sound box-Star of the attention weights A = softmax(logits) where
            % S_logits is a matrix-Star (dim = prod(shape)) of pre-softmax logits,
            % shape = [R n], softmax over the n key axis. Concretizes logits to a
            % sound box then applies the exact correlated row bound. Predicates are
            % NOT preserved (this is the IBP-class weight set the value path consumes).
            if nargin < 3 || isempty(mode), mode = 'estimate'; end
            R = shape(1); n = shape(2);
            if S_logits.dim ~= R*n, error('SoftmaxAttn:shape','logits dim mismatch'); end
            [lo, hi] = SoftmaxAttn.starBounds(S_logits, mode);
            s_lo = reshape(lo, [R n]);
            s_hi = reshape(hi, [R n]);
            [a_lb, a_ub] = SoftmaxAttn.correlatedRowSoftmaxBounds(s_lo, s_hi);
            W = SoftmaxAttn.boxStar(reshape(a_lb, [], 1), reshape(a_ub, [], 1));
        end

        % ----- symbolic A*V envelope (value-path-preserving) ---------------
        function O = avEnvelopeStar(a_lb, a_ub, V_star, shV, mode)
            % Sound symbolic reach of O = A * V where A in [a_lb,a_ub] (the softmax
            % weights, a_lb>=0) is a fixed interval and V is a SYMBOLIC Star whose
            % column-major reshape is [K x D] (shV=[K D]). a_lb,a_ub are [M x K].
            % Returns a Star O of dim M*D, AFFINE IN V's PREDICATES (so O stays
            % correlated with the network input). V's predicates are kept as a
            % prefix; M*D fresh predicates (one per output) are appended, each
            % sandwiched by two sign-aware McCormick facets.
            %
            % Per term A[m,k]*V[k,d] with a in [al,au] (>=0), v in [vlo,vhi]:
            %   v>=0  : up = au*v,             low = al*v
            %   v<=0  : up = al*v,             low = au*v      (multipliers swap)
            %   mixed : up = chord of max_a(a*v) over [vlo,vhi]  (convex),
            %           low = chord of min_a(a*v) over [vlo,vhi] (concave)
            % All slopes are >=0 because a_lb,a_ub>=0, so the per-output box is
            % o_ub = sum_k up_slope*vhi + up_bias,  o_lb = sum_k low_slope*vlo + low_bias.
            if nargin < 5 || isempty(mode), mode = 'estimate'; end
            if any(a_lb(:) < -1e-12)
                error('SoftmaxAttn:avNegativeWeight', ...
                    'avEnvelopeStar requires a_lb>=0 (softmax weights); got min=%g', min(a_lb(:)));
            end
            K = shV(1); D = shV(2);
            [M, Ka] = size(a_lb);
            if Ka ~= K, error('SoftmaxAttn:shape','A inner dim %d ~= K %d', Ka, K); end
            if V_star.dim ~= K*D, error('SoftmaxAttn:shape','V dim %d ~= K*D %d', V_star.dim, K*D); end

            n_old = V_star.nVar;
            cV = reshape(V_star.V(:,1), [K, D]);
            if n_old > 0
                GV = reshape(V_star.V(:,2:end), [K, D, n_old]);   % generators
            else
                GV = zeros(K, D, 0);
            end
            [vlo_v, vhi_v] = SoftmaxAttn.starBounds(V_star, mode);
            vlo = reshape(vlo_v, [K, D]);
            vhi = reshape(vhi_v, [K, D]);

            % sign-aware envelope slopes/biases as [M,K,D] via implicit expansion
            AL = repmat(a_lb, [1, 1, D]);                 % [M,K,D]
            AU = repmat(a_ub, [1, 1, D]);
            VLO = repmat(reshape(vlo, [1, K, D]), [M, 1, 1]);
            VHI = repmat(reshape(vhi, [1, K, D]), [M, 1, 1]);
            pos = VLO >= 0;
            neg = VHI <= 0;
            mixed = ~pos & ~neg;

            up_slope = zeros(M,K,D); up_bias = zeros(M,K,D);
            low_slope = zeros(M,K,D); low_bias = zeros(M,K,D);
            % pos branch
            up_slope(pos) = AU(pos);   low_slope(pos) = AL(pos);
            % neg branch (swap)
            up_slope(neg) = AL(neg);   low_slope(neg) = AU(neg);
            % mixed branch (McCormick secant)
            w = VHI - VLO;  w(w < 1e-12) = 1e-12;
            us = (AU.*VHI - AL.*VLO) ./ w;  ub = AL.*VLO - us.*VLO;
            ls = (AL.*VHI - AU.*VLO) ./ w;  lb = AU.*VLO - ls.*VLO;
            up_slope(mixed) = us(mixed);   up_bias(mixed) = ub(mixed);
            low_slope(mixed) = ls(mixed);  low_bias(mixed) = lb(mixed);
            % FP guard: slopes are >=0 by construction (a>=0)
            up_slope = max(up_slope, 0);  low_slope = max(low_slope, 0);

            % contract over k -> coef[m,d,:], const[m,d], box o_lb/o_ub[m,d]
            n_out = M*D;
            CoefUp  = zeros(n_out, n_old);
            CoefLow = zeros(n_out, n_old);
            const_up  = zeros(M, D);
            const_low = zeros(M, D);
            o_ub = zeros(M, D);
            o_lb = zeros(M, D);
            for m = 1:M
                for d = 1:D
                    us_ = reshape(up_slope(m,:,d),  [K,1]);
                    ls_ = reshape(low_slope(m,:,d), [K,1]);
                    ub_ = reshape(up_bias(m,:,d),   [K,1]);
                    lb_ = reshape(low_bias(m,:,d),  [K,1]);
                    Gkd = reshape(GV(:,d,:), [K, n_old]);   % generators of V[:,d]
                    idx = (d-1)*M + m;                      % column-major O[m,d]
                    if n_old > 0
                        CoefUp(idx, :)  = us_' * Gkd;
                        CoefLow(idx, :) = ls_' * Gkd;
                    end
                    cVd = cV(:, d);
                    const_up(m,d)  = us_' * cVd + sum(ub_);
                    const_low(m,d) = ls_' * cVd + sum(lb_);
                    o_ub(m,d) = us_' * vhi(:,d) + sum(ub_);
                    o_lb(m,d) = ls_' * vlo(:,d) + sum(lb_);
                end
            end
            const_up_v  = reshape(const_up,  [n_out, 1]);
            const_low_v = reshape(const_low, [n_out, 1]);
            o_ub_v = reshape(o_ub, [n_out, 1]);
            o_lb_v = reshape(o_lb, [n_out, 1]);

            % assemble output Star: state O = fresh predicates (center 0)
            Vout = [zeros(n_out,1), zeros(n_out, n_old), eye(n_out)];
            Cold = [V_star.C, zeros(size(V_star.C,1), n_out)];
            Cup  = [-CoefUp,  eye(n_out)];     %  o - coef_up*a  <=  const_up
            Clow = [ CoefLow, -eye(n_out)];    % coef_low*a - o  <= -const_low
            C = [Cold; Cup; Clow];
            dvec = [V_star.d; const_up_v; -const_low_v];
            plb = [V_star.predicate_lb; o_lb_v];
            pub = [V_star.predicate_ub; o_ub_v];
            O = Star(Vout, C, dvec, plb, pub);
        end

        % ----- single-head attention kernel --------------------------------
        function O = singleHeadAttn(Q, K, V, scale, shape, mode)
            % Sound reach of O = softmax(scale * Q*K') * V for one head.
            % Q,K are matrix-Stars of dim N*D (column-major reshape [N D]) used only
            % for sound score bounds. V is a matrix-Star of dim N*Dv (its value
            % width Dv = d_v may differ from the query/key width D = d_k; Dv is
            % inferred from V). V is kept SYMBOLIC, so O (dim N*Dv) stays correlated
            % with V's predicates (value path). softmax over keys.
            if nargin < 6 || isempty(mode), mode = 'estimate'; end
            N = shape(1); D = shape(2);
            [ql, qu] = SoftmaxAttn.starBounds(Q, mode);
            [kl, ku] = SoftmaxAttn.starBounds(K, mode);
            Ql = reshape(ql,[N D]); Qu = reshape(qu,[N D]);
            Kl = reshape(kl,[N D]); Ku = reshape(ku,[N D]);
            % scores = Q * K' : intervalMatMul with the second operand transposed
            [scl, scu] = SoftmaxAttn.intervalMatMul(Ql, Qu, Kl', Ku');   % [N N]
            lo = min(scale.*scl, scale.*scu);
            hi = max(scale.*scl, scale.*scu);
            [a_lb, a_ub] = SoftmaxAttn.correlatedRowSoftmaxBounds(lo, hi);
            % The value width may differ from the query/key width (d_v ~= d_k); the
            % score path used D above, but the A*V envelope must shape V by its own
            % width, inferred from V (= d_v). For the ViT case d_v == d_k == D.
            Dv = round(V.dim / N);
            if Dv * N ~= V.dim
                error('SoftmaxAttn:shape', 'V dim %d not a multiple of N=%d', V.dim, N);
            end
            O = SoftmaxAttn.avEnvelopeStar(a_lb, a_ub, V, [N Dv], mode);
        end

        % ----- per-token linear map (exact, predicate-preserving) ----------
        function Y = perTokenLinear(X, M, b, N)
            % Apply Y = Xt*M + b per token, where X is a matrix-Star of dim N*Ein
            % (e-major: reshape(state,[N Ein])) and M is [Ein x Eout], b is [Eout x 1].
            % vec(Xt*M) = kron(M', I_N) vec(Xt); bias broadcast over tokens.
            Ein = size(M, 1); Eout = size(M, 2);
            if X.dim ~= N*Ein, error('SoftmaxAttn:perTokenLinear','dim %d ~= N*Ein %d', X.dim, N*Ein); end
            W = kron(M', eye(N));
            bb = kron(b(:), ones(N, 1));
            Y = X.affineMap(W, bb);
        end

        % ----- slice one head's channels from a matrix-Star -----------------
        function Sh = sliceHead(S, h, D, N)
            % Channels (h-1)*D+1 : h*D of an e-major matrix-Star of dim N*E.
            idx = ((h-1)*D)*N + (1:D*N);
            Sh = Star(S.V(idx, :), S.C, S.d, S.predicate_lb, S.predicate_ub);
        end

        % ----- full multi-head self-attention reach (no residual) -----------
        function AO = selfAttentionReach(X, N, P, mode)
            % Sound reach of one ViT encoder attention sublayer:
            %   Q,K,V = perTokenLinear(X);  per head softmax(scale*Qh*Kh')*Vh;
            %   concat heads;  out = perTokenLinear(O).
            % X: e-major matrix-Star of dim N*E. P holds E,H,D,scale and the four
            % projection weights/biases (Mq,bq,Mk,bk,Mv,bv,Mo,bo) as ONNX MatMul
            % matrices ([Ein x Eout]). Returns AO (dim N*E), V-path symbolic.
            if nargin < 4 || isempty(mode), mode = 'estimate'; end
            E = P.E; H = P.H; D = P.D; scale = P.scale;
            n_old = X.nVar;
            Q = SoftmaxAttn.perTokenLinear(X, P.Mq, P.bq, N);
            K = SoftmaxAttn.perTokenLinear(X, P.Mk, P.bk, N);
            V = SoftmaxAttn.perTokenLinear(X, P.Mv, P.bv, N);
            Os = cell(1, H);
            for h = 1:H
                Qh = SoftmaxAttn.sliceHead(Q, h, D, N);
                Kh = SoftmaxAttn.sliceHead(K, h, D, N);
                Vh = SoftmaxAttn.sliceHead(V, h, D, N);
                Os{h} = SoftmaxAttn.singleHeadAttn(Qh, Kh, Vh, scale, [N D], mode);
            end
            O = SoftmaxAttn.prefixConcat(Os, n_old);   % e-major [N E] layout
            AO = SoftmaxAttn.perTokenLinear(O, P.Mo, P.bo, N);
        end

        % ----- prefix-aligned residual add ---------------------------------
        function S = prefixAdd(x, branch)
            % Sound, tight residual y = x + branch where x's predicates are a
            % genuine PREFIX of branch's (branch was built from x by appending
            % relaxation predicates). Adds over branch's (superset) constraint
            % system, so x's predicate correlation is preserved exactly -- no
            % block-diagonal Minkowski blow-up, and never the structural-equality
            % "shared alpha" fast path that under-approximates independent operands.
            if x.dim ~= branch.dim
                error('SoftmaxAttn:prefixAdd', 'dim mismatch (%d vs %d)', x.dim, branch.dim);
            end
            nx = x.nVar; nb = branch.nVar;
            if nx > nb
                error('SoftmaxAttn:prefixAdd', 'x (nVar=%d) is not a prefix of branch (nVar=%d)', nx, nb);
            end
            Vx = zeros(x.dim, nb + 1);
            Vx(:, 1:nx+1) = x.V;             % center + x's prefix generators
            Vnew = Vx + branch.V;
            S = Star(Vnew, branch.C, branch.d, branch.predicate_lb, branch.predicate_ub);
        end

        % ----- prefix-aligned concatenation (assemble heads) ---------------
        function S = prefixConcat(stars, n_old)
            % Concatenate the STATES of several stars that share the same first
            % n_old predicates (the input prefix) into one star whose fresh tails
            % are block-diagonal. Sound: the result is the intersection of each
            % star's constraint system (all sharing the prefix), with disjoint
            % fresh predicates. Used to assemble attention heads / Q,K,V streams.
            H = numel(stars);
            tails = zeros(1, H); dims = zeros(1, H);
            for h = 1:H
                tails(h) = stars{h}.nVar - n_old;
                dims(h)  = stars{h}.dim;
                if tails(h) < 0
                    error('SoftmaxAttn:prefixConcat', 'star %d has < n_old predicates', h);
                end
            end
            total_dim  = sum(dims);
            total_tail = sum(tails);
            total_pred = n_old + total_tail;

            Vnew = zeros(total_dim, total_pred + 1);
            Crows = {}; drows = {};
            tail_plb = []; tail_pub = [];
            plb_prefix = stars{1}.predicate_lb(1:n_old);
            pub_prefix = stars{1}.predicate_ub(1:n_old);

            rofs = 0; tofs = 0;
            for h = 1:H
                Sh = stars{h};
                dh = dims(h); th = tails(h);
                rows = rofs + (1:dh);
                Vnew(rows, 1) = Sh.V(:, 1);                       % center
                if n_old > 0
                    Vnew(rows, 2:n_old+1) = Sh.V(:, 2:n_old+1);   % shared-prefix generators
                end
                if th > 0
                    Vnew(rows, n_old+1+tofs + (1:th)) = Sh.V(:, n_old+1 + (1:th));  % tail block
                    tail_plb = [tail_plb; Sh.predicate_lb(n_old + (1:th))];
                    tail_pub = [tail_pub; Sh.predicate_ub(n_old + (1:th))];
                end
                Ch = Sh.C;
                if ~isempty(Ch)
                    nc = size(Ch, 1);
                    Cmap = zeros(nc, total_pred);
                    Cmap(:, 1:n_old) = Ch(:, 1:n_old);
                    if th > 0
                        Cmap(:, n_old+tofs + (1:th)) = Ch(:, n_old + (1:th));
                    end
                    if h == 1
                        Crows{end+1} = Cmap;  drows{end+1} = Sh.d;   % keep shared-prefix rows once
                    else
                        if th > 0
                            keep = any(Ch(:, n_old + (1:th)) ~= 0, 2);  % only this head's tail facets
                        else
                            keep = false(nc, 1);
                        end
                        Crows{end+1} = Cmap(keep, :);  drows{end+1} = Sh.d(keep);
                    end
                end
                rofs = rofs + dh; tofs = tofs + th;
            end
            Cnew = vertcat(Crows{:});
            dnew = vertcat(drows{:});
            plb = [plb_prefix; tail_plb];
            pub = [pub_prefix; tail_pub];
            S = Star(Vnew, Cnew, dnew, plb, pub);
        end

        % ----- robust box-Star constructor ---------------------------------
        function S = boxStar(lo, hi)
            % Build a box Star from [lo,hi] column vectors, tolerant of
            % degenerate (lo==hi) coordinates. Mirrors n2v Star.from_bounds:
            % center c, one generator per coordinate (radius r), explicit
            % predicate box alpha in [-1,1] and the +/- I constraint rows.
            lo = double(lo(:)); hi = double(hi(:));
            d = numel(lo);
            % numerical guard: clamp tiny inversions from FP
            sw = hi < lo;
            if any(sw), tmp = lo(sw); lo(sw) = hi(sw); hi(sw) = tmp; end
            c = (lo + hi) / 2;
            r = (hi - lo) / 2;
            V = [c, diag(r)];
            C = [eye(d); -eye(d)];
            dvec = ones(2*d, 1);
            plb = -ones(d,1); pub = ones(d,1);
            S = Star(V, C, dvec, plb, pub);
        end

    end
end
