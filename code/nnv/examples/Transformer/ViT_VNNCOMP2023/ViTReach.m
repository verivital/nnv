classdef ViTReach
    % ViTReach - sound star-set reachability driver for the VNN-COMP 2023 ViT
    % benchmark (CIFAR-10, L-inf eps=1/255, argmax-preservation), implemented
    % entirely on NNV Star sets via the SoftmaxAttn sound primitives.
    %
    % The benchmark ONNX (opset 9, Slice v1 + per-token BatchNorm) is not ingested
    % by MATLAB's importNetworkFromONNX, so -- as in the n2v reference -- we mirror
    % the exported graph explicitly and load weights from a .mat produced by
    % extract_weights.py. evaluate() is parity-checked to ~1e-5 vs onnxruntime;
    % reach() is the sound over-approximation that propagates a per-pixel eps-box.
    %
    % Architecture (ViT_BN): Conv patch-embed -> [cls] -> +pos ->
    %   depth x ( x = x + Attn(BN(x));  x = x + FF_ReLU(BN(x)) ) ->
    %   ReduceMean over tokens -> BN -> Linear head (10 logits).
    % Only genuine nonlinearities: softmax attention (bilinear QK^T + softmax + A*V)
    % and the FF ReLU. Everything else folds to exact affine maps.
    %
    % Token-state convention: a matrix-Star of dim N*E in e-major layout, i.e.
    % reshape(state,[N E]) recovers the [token x channel] matrix. See SoftmaxAttn.
    %
    % Author: NNV Team (ViT track). Date 2026.

    methods (Static)

        % ================= model loading =================
        function M = load(matfile)
            S = load(matfile);
            M = S;
            M.patch = double(S.patch); M.depth = double(S.depth);
            M.dim = double(S.dim); M.heads = double(S.heads); M.dim_head = double(S.dim_head);
            M.scale = 1 / sqrt(M.dim_head);
            M.Hg = 32 / M.patch;
            M.n_patch = M.Hg^2;
            M.N = M.n_patch + 1;
            M.cls = double(S.cls(:));
            M.pos = double(S.pos);              % [N x E]
            M.proj_w = double(S.proj_w);        % [E,3,p,p]
            M.proj_b = double(S.proj_b(:));
            M.head_W = double(S.head_W);        % [10 x E]
            M.head_b = double(S.head_b(:));
            M.head_bn_scale = double(S.head_bn_scale(:));
            M.head_bn_shift = double(S.head_bn_shift(:));
            % blocks come in as a cell array of structs (scipy savemat of a list);
            % normalise field orientations into a clean struct array.
            Bc = S.blocks;
            if ~iscell(Bc), Bc = num2cell(Bc); end
            for i = 1:numel(Bc)
                b = Bc{i};
                bb(i).bn1_scale = double(b.bn1_scale(:)); bb(i).bn1_shift = double(b.bn1_shift(:));
                bb(i).bn2_scale = double(b.bn2_scale(:)); bb(i).bn2_shift = double(b.bn2_shift(:));
                bb(i).Mq = double(b.Mq); bb(i).bq = double(b.bq(:));
                bb(i).Mk = double(b.Mk); bb(i).bk = double(b.bk(:));
                bb(i).Mv = double(b.Mv); bb(i).bv = double(b.bv(:));
                bb(i).Mo = double(b.Mo); bb(i).bo = double(b.bo(:));
                bb(i).ff1_W = double(b.ff1_W); bb(i).ff1_b = double(b.ff1_b(:));
                bb(i).ff2_W = double(b.ff2_W); bb(i).ff2_b = double(b.ff2_b(:));
            end
            M.blocks = bb;
        end

        % ================= concrete forward (parity oracle) =================
        function logits = evaluate(M, img)
            % img: [3,32,32] already normalised. Returns logits [10x1].
            p = M.patch; Hg = M.Hg; E = M.dim; H = M.heads; D = M.dim_head;
            z = zeros(M.n_patch, E);
            for gh = 0:Hg-1
                for gw = 0:Hg-1
                    p0 = gh*Hg + gw;                                  % row-major patch index
                    patch = img(:, gh*p+(1:p), gw*p+(1:p));           % [3,p,p]
                    for e = 1:E
                        w = reshape(M.proj_w(e,:,:,:), [3 p p]);
                        z(p0+1, e) = M.proj_b(e) + sum(w(:).*patch(:));
                    end
                end
            end
            Z = [M.cls'; z] + M.pos;                                  % [N x E]
            for i = 1:M.depth
                blk = M.blocks(i);
                Z = Z + ViTReach.attnFwd(ViTReach.bnFwd(Z, blk.bn1_scale, blk.bn1_shift), blk, H, D, M.scale);
                Z = Z + ViTReach.ffFwd(ViTReach.bnFwd(Z, blk.bn2_scale, blk.bn2_shift), blk);
            end
            zt = mean(Z, 1)';                                         % [E x 1]
            zt = zt .* M.head_bn_scale + M.head_bn_shift;
            logits = M.head_W * zt + M.head_b;                        % [10 x 1]
        end

        function Y = bnFwd(Z, scale, shift)
            Y = Z .* scale' + shift';                                 % [N x E] per-channel affine
        end

        function out = attnFwd(X, blk, H, D, scale)
            [N, E] = size(X);
            Q = X*blk.Mq + blk.bq';  K = X*blk.Mk + blk.bk';  V = X*blk.Mv + blk.bv';
            O = zeros(N, E);
            for h = 1:H
                cols = (h-1)*D + (1:D);
                A = ViTReach.rowSoftmax(scale*(Q(:,cols)*K(:,cols)'));
                O(:, cols) = A * V(:, cols);
            end
            out = O*blk.Mo + blk.bo';
        end

        function out = ffFwd(X, blk)
            out = max(X*blk.ff1_W + blk.ff1_b', 0) * blk.ff2_W + blk.ff2_b';
        end

        function A = rowSoftmax(S)
            m = max(S,[],2); ex = exp(S-m); A = ex./sum(ex,2);
        end

        % ================= sound reachability =================
        function L = reach(M, lb, ub, opt)
            % lb,ub: per-pixel normalised bounds, ordered (c,h,w) with
            % idx=(c-1)*1024+(h-1)*32+(w-1)+1. Returns logits Star (dim 10).
            % opt fields: .mode ('estimate'|'lp', default 'estimate'),
            %             .relu ('approx-star' default).
            if nargin < 4, opt = struct(); end
            if ~isfield(opt,'mode'), opt.mode = 'estimate'; end
            if ~isfield(opt,'relu'), opt.relu = 'fast'; end
            E = M.dim; H = M.heads; D = M.dim_head; N = M.N;

            Xin = Star(lb(:), ub(:));                                 % dim 3072
            [Wpe, bpe] = ViTReach.patchEmbed(M);
            X = Xin.affineMap(Wpe, bpe);                              % token-state dim N*E
            for i = 1:M.depth
                blk = M.blocks(i);
                % attention sublayer
                P = ViTReach.attnParams(blk, E, H, M.scale);
                Xn = ViTReach.bnStar(X, blk.bn1_scale, blk.bn1_shift, N);
                AO = SoftmaxAttn.selfAttentionReach(Xn, N, P, opt.mode);
                X = SoftmaxAttn.prefixAdd(X, AO);
                % feed-forward sublayer
                Xn = ViTReach.bnStar(X, blk.bn2_scale, blk.bn2_shift, N);
                F1 = SoftmaxAttn.perTokenLinear(Xn, blk.ff1_W, blk.ff1_b, N);   % dim N*96
                R  = ViTReach.reluStar(F1, opt.relu);
                F2 = SoftmaxAttn.perTokenLinear(R, blk.ff2_W, blk.ff2_b, N);    % dim N*E
                X = SoftmaxAttn.prefixAdd(X, F2);
            end
            % ReduceMean over tokens -> dim E
            Wm = ViTReach.meanMap(N, E);
            zt = X.affineMap(Wm, zeros(E,1));
            % head BN (per-channel affine) + head linear
            zt = zt.affineMap(diag(M.head_bn_scale), M.head_bn_shift);
            L = zt.affineMap(M.head_W, M.head_b);                     % dim 10
        end

        function [robust, margins, L] = verify(M, lb, ub, label, opt)
            % Returns robust in {1 robust, 2 unknown}. margins(i) = lower bound on
            % (logit_label - logit_i) for i~=label; robust iff all > 0.
            if nargin < 5, opt = struct(); end
            if ~isfield(opt,'marginMode'), opt.marginMode = 'estimate'; end
            L = ViTReach.reach(M, lb, ub, opt);
            label = label + 1;                                        % 0-based -> 1-based
            margins = inf(10,1);
            robust = 1;
            for i = 1:10
                if i == label, continue; end
                c = zeros(10,1); c(label) = 1; c(i) = -1;
                ms = L.affineMap(c', 0);                              % scalar Star = logit_L - logit_i
                if strcmpi(opt.marginMode,'lp')
                    mlb = ms.getMin(1, 'linprog');
                else
                    [mlb, ~] = ms.estimateRanges();
                end
                margins(i) = mlb;
                if mlb <= 0, robust = 2; end
            end
        end

        % ---- reach helpers ----
        function Y = bnStar(X, scale, shift, N)
            d = kron(scale(:), ones(N,1));
            b = kron(shift(:), ones(N,1));
            Y = X.affineMap(diag(d), b);
        end

        function R = reluStar(S, method)
            % Sound ReLU on a Star.
            %   'fast' (default): symbolic triangle relaxation using cheap ESTIMATE
            %     neuron bounds -- no LP, preserves the predicate prefix (appends one
            %     fresh predicate per unstable neuron). This is what makes the forward
            %     reach run in seconds rather than minutes (LP neuron bounds over the
            %     few-thousand-predicate residual stream dominate otherwise).
            %   otherwise: delegate to NNV's PosLin func (LP-tight but slow).
            if nargin < 2 || isempty(method) || strcmpi(method, 'fast')
                R = ViTReach.reluStarApproxFast(S);
            else
                Rs = PosLin.reach(S, method);
                R = Rs(1);
            end
        end

        function R = reluStarApproxFast(S)
            % approx-star ReLU with estimate (box) neuron bounds; sound triangle.
            [l, u] = S.estimateRanges();
            n = S.dim; nVar = S.nVar;
            c = S.V(:, 1); G = S.V(:, 2:end);
            unstable = find(l < 0 & u > 0);
            posn = find(l >= 0);
            nr = numel(unstable);

            newV = zeros(n, nVar + 1 + nr);
            newV(posn, 1) = c(posn);
            if nVar > 0, newV(posn, 2:nVar+1) = G(posn, :); end   % positive-stable: identity
            % negative-stable rows stay 0 (ReLU outputs 0)
            for k = 1:nr
                newV(unstable(k), nVar+1+k) = 1;                  % unstable: y = fresh pred
            end

            if isempty(S.C)
                Cold = zeros(0, nVar + nr); dold = zeros(0, 1);
            else
                Cold = [S.C, zeros(size(S.C,1), nr)]; dold = S.d;
            end
            plb = [S.predicate_lb; zeros(nr,1)];
            pub = [S.predicate_ub; zeros(nr,1)];

            rows = zeros(3*nr, nVar + nr); rd = zeros(3*nr, 1);
            for k = 1:nr
                i = unstable(k); li = l(i); ui = u(i);
                xc = c(i); xg = G(i, :);
                col = nVar + k;
                lam = ui / (ui - li);
                % a >= 0
                rows(3*k-2, col) = -1;                       rd(3*k-2) = 0;
                % a >= x_i   (x_i - a <= 0)
                rows(3*k-1, 1:nVar) = xg; rows(3*k-1, col) = -1;  rd(3*k-1) = -xc;
                % a <= lam*(x_i - li)  ( -lam*xg*alpha + a <= lam*(xc - li) )
                rows(3*k,   1:nVar) = -lam*xg; rows(3*k, col) = 1; rd(3*k) = lam*(xc - li);
                pub(col) = ui;
            end
            R = Star(newV, [Cold; rows], [dold; rd], plb, pub);
        end

        function Wm = meanMap(N, E)
            % [E x N*E] averaging over tokens of an e-major matrix-Star
            Wm = zeros(E, N*E);
            for e = 1:E
                Wm(e, (e-1)*N + (1:N)) = 1/N;
            end
        end

        function P = attnParams(blk, E, H, scale)
            P.E=E; P.H=H; P.D=E/H; P.scale=scale;
            P.Mq=blk.Mq; P.bq=blk.bq; P.Mk=blk.Mk; P.bk=blk.bk;
            P.Mv=blk.Mv; P.bv=blk.bv; P.Mo=blk.Mo; P.bo=blk.bo;
        end

        function [Wpe, bpe] = patchEmbed(M)
            % exact affine map from the (c,h,w)-ordered pixel vector (dim 3072) to
            % the e-major token-state (dim N*E), folding conv patch-embed + cls + pos.
            p = M.patch; Hg = M.Hg; E = M.dim; N = M.N;
            Wpe = zeros(N*E, 3*32*32);
            bpe = zeros(N*E, 1);
            for e = 1:E
                bpe((e-1)*N + 1) = M.cls(e) + M.pos(1, e);           % cls token (token 1)
            end
            for gh = 0:Hg-1
                for gw = 0:Hg-1
                    p0 = gh*Hg + gw;
                    n = p0 + 2;                                       % token row (cls is 1)
                    for e = 1:E
                        bpe((e-1)*N + n) = M.proj_b(e) + M.pos(n, e);
                        for c = 1:3
                            for ii = 0:p-1
                                for jj = 0:p-1
                                    h1 = gh*p + ii + 1; w1 = gw*p + jj + 1;
                                    col = (c-1)*1024 + (h1-1)*32 + (w1-1) + 1;
                                    Wpe((e-1)*N + n, col) = M.proj_w(e, c, ii+1, jj+1);
                                end
                            end
                        end
                    end
                end
            end
        end

        % ================= input box construction =================
        function [lb, ub] = epsBox(M, img01, eps)
            % img01: [3,32,32] raw pixels in [0,1]. Returns normalised per-pixel
            % bounds in the (c,h,w) ordering used by patchEmbed.
            if nargin < 3, eps = 1/255; end
            mean_ = M.mean(:); std_ = M.std(:);
            lo = min(max(img01 - eps, 0), 1);
            hi = min(max(img01 + eps, 0), 1);
            lb = zeros(3*32*32,1); ub = zeros(3*32*32,1);
            for c = 1:3
                for h = 1:32
                    for w = 1:32
                        idx = (c-1)*1024 + (h-1)*32 + (w-1) + 1;
                        lb(idx) = (lo(c,h,w) - mean_(c)) / std_(c);
                        ub(idx) = (hi(c,h,w) - mean_(c)) / std_(c);
                    end
                end
            end
        end

        function x = normImg(M, img01)
            % normalise a raw [3,32,32] image to the model's input space
            x = (img01 - reshape(M.mean,[3 1 1])) ./ reshape(M.std,[3 1 1]);
        end

    end
end
