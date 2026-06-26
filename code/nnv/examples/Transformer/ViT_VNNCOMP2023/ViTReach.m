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
        function [L, info] = reach(M, lb, ub, opt)
            % lb,ub: per-pixel normalised bounds, ordered (c,h,w) with
            % idx=(c-1)*1024+(h-1)*32+(w-1)+1. Returns logits Star (dim 10).
            % opt fields: .mode ('estimate'|'lp', default 'estimate'),
            %             .relu ('fast' default),
            %             .splits (struct array .block/.neuron/.phase forcing FF
            %                      ReLU phases for branch-and-bound; +1 active, -1 inactive).
            % info.reluBounds{i} = [l u] of block i's FF ReLU input (for BaB choice).
            if nargin < 4, opt = struct(); end
            if ~isfield(opt,'mode'), opt.mode = 'estimate'; end
            if ~isfield(opt,'relu'), opt.relu = 'fast'; end
            if ~isfield(opt,'splits'), opt.splits = []; end
            E = M.dim; H = M.heads; D = M.dim_head; N = M.N;
            info = struct('reluBounds', {cell(1, M.depth)});

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
                [lF, uF] = F1.estimateRanges(); info.reluBounds{i} = [lF(:), uF(:)];
                forced = ViTReach.forcedFor(opt.splits, i);
                R  = ViTReach.reluStar(F1, opt.relu, forced);
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

        function forced = forcedFor(splits, blockIdx)
            forced = zeros(0,2);
            if isempty(splits), return; end
            for s = 1:numel(splits)
                if splits(s).block == blockIdx
                    forced(end+1,:) = [splits(s).neuron, splits(s).phase]; %#ok<AGROW>
                end
            end
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

        % ================= branch-and-bound =================
        function [status, info] = verifyBaB(M, img, label, eps, opt)
            % Sound input-splitting branch-and-bound robustness verification.
            %   status: 1 robust, 0 not-robust (counterexample), 2 unknown (budget).
            % Uses ViTReach.reach as the bounding oracle (the new sound ViT layers)
            % and bisects the per-pixel eps-box on the widest dimension. Sound:
            % a box is robust only if every leaf is certified; any falsifying
            % sample in any sub-box proves not-robust. Complete in the limit
            % (subject to opt.babMaxNodes). Input BaB is weak in 3072-d (the
            % effective gap-closer is neuron/ReLU-split beta-CROWN; see STATUS.md),
            % but this demonstrates a sound BaB driven by the attention layers and
            % can certify a strictly larger radius than single-shot reach.
            if nargin < 5, opt = struct(); end
            if ~isfield(opt,'mode'), opt.mode = 'estimate'; end
            if ~isfield(opt,'relu'), opt.relu = 'fast'; end
            if ~isfield(opt,'marginMode'), opt.marginMode = 'estimate'; end
            if ~isfield(opt,'babMaxNodes'), opt.babMaxNodes = 50; end
            if ~isfield(opt,'babFalsifyTries'), opt.babFalsifyTries = 100; end

            [lb0, ub0] = ViTReach.epsBox(M, img, eps);
            stack = {struct('lb',lb0,'ub',ub0)};
            nodes = 0; anyUnknown = false;
            while ~isempty(stack)
                node = stack{end}; stack(end) = [];
                nodes = nodes + 1;
                [rob, margins] = ViTReach.verifyBox(M, node.lb, node.ub, label, opt);
                if rob == 1
                    continue;                                  % leaf certified robust
                end
                % try to falsify this sub-box (sound not-robust witness)
                ce = ViTReach.falsifyBox(M, node.lb, node.ub, label, opt.babFalsifyTries);
                if ce
                    status = 0; info = struct('nodes',nodes,'reason','counterexample');
                    return;
                end
                if nodes >= opt.babMaxNodes
                    anyUnknown = true; continue;               % budget hit: leaf unknown
                end
                % split the widest dimension into two halves
                [~, d] = max(node.ub - node.lb);
                mid = (node.lb(d) + node.ub(d)) / 2;
                left = node;  left.ub(d) = mid;
                right = node; right.lb(d) = mid;
                stack{end+1} = left; stack{end+1} = right; %#ok<AGROW>
            end
            if anyUnknown, status = 2; else, status = 1; end
            info = struct('nodes',nodes,'reason','exhausted');
        end

        function [status, info] = verifyBaBRelu(M, img, label, eps, opt)
            % Sound FF-ReLU phase-split branch-and-bound (beta-CROWN style) driven
            % by the sound ViT layers. A node forces a set of ReLU phases; its two
            % children fix one more candidate neuron active/inactive (their union is
            % the node), so the node is robust iff both children are. The binding
            % margin is tightened with an LP that honours the forced half-spaces.
            %   status: 1 robust, 2 unknown (budget). (No CE search here -- the input
            %   box is fixed; pair with verifyBaB/falsify for not-robust witnesses.)
            if nargin < 5, opt = struct(); end
            if ~isfield(opt,'mode'), opt.mode = 'estimate'; end
            if ~isfield(opt,'relu'), opt.relu = 'fast'; end
            if ~isfield(opt,'babMaxNodes'), opt.babMaxNodes = 40; end
            if ~isfield(opt,'maxCand'), opt.maxCand = 14; end
            [lb, ub] = ViTReach.epsBox(M, img, eps);

            % root reach: get candidate unstable FF-ReLU neurons + instability scores
            opt0 = opt; opt0.splits = [];
            [~, info0] = ViTReach.reach(M, lb, ub, opt0);
            cand = zeros(0,3);   % [block neuron score]
            for i = 1:M.depth
                lu = info0.reluBounds{i};
                uns = find(lu(:,1) < 0 & lu(:,2) > 0);
                sc = min(-lu(uns,1), lu(uns,2));
                cand = [cand; [i*ones(numel(uns),1), uns, sc]]; %#ok<AGROW>
            end
            [~, ord] = sort(cand(:,3), 'descend');
            cand = cand(ord, :);
            T = min(opt.maxCand, size(cand,1));
            cand = cand(1:T, :);

            stack = {struct('splits', struct('block',{},'neuron',{},'phase',{}), 'depth', 0)};
            nodes = 0; anyUnknown = false;
            while ~isempty(stack)
                node = stack{end}; stack(end) = [];
                nodes = nodes + 1;
                rob = ViTReach.verifyBoxSplits(M, lb, ub, label, opt, node.splits);
                if rob == 1, continue; end                 % leaf certified
                if node.depth >= T || nodes >= opt.babMaxNodes
                    anyUnknown = true; continue;           % budget / no more candidates
                end
                c = cand(node.depth+1, :);
                for ph = [1, -1]
                    s = node.splits;
                    s(end+1) = struct('block',c(1),'neuron',c(2),'phase',ph); %#ok<AGROW>
                    stack{end+1} = struct('splits', s, 'depth', node.depth+1); %#ok<AGROW>
                end
            end
            if anyUnknown, status = 2; else, status = 1; end
            info = struct('nodes', nodes, 'nCand', T);
        end

        function rob = verifyBoxSplits(M, lb, ub, label, opt, splits)
            % robustness of [lb,ub] under forced ReLU phases: estimate margins
            % first (cheap), then an LP (honouring the split half-spaces) on any
            % class the estimate cannot clear. rob = 1 robust, 2 unknown.
            o = opt; o.splits = splits;
            L = ViTReach.reach(M, lb, ub, o);
            lab = label + 1; rob = 1;
            for i = 1:10
                if i == lab, continue; end
                cc = zeros(10,1); cc(lab) = 1; cc(i) = -1;
                ms = L.affineMap(cc', 0);
                [e,~] = ms.estimateRanges();
                if e > 0, continue; end                    % cleared by estimate
                if ms.getMin(1,'linprog') > 0, continue; end  % cleared by LP (uses splits)
                rob = 2; return;
            end
        end

        function [rob, margins] = verifyBox(M, lb, ub, label, opt)
            % verify robustness on an explicit normalised pixel box [lb,ub].
            L = ViTReach.reach(M, lb, ub, opt);
            lab = label + 1;
            margins = inf(10,1); rob = 1;
            for i = 1:10
                if i == lab, continue; end
                c = zeros(10,1); c(lab) = 1; c(i) = -1;
                ms = L.affineMap(c', 0);
                if strcmpi(opt.marginMode,'lp'), mlb = ms.getMin(1,'linprog');
                else, [mlb,~] = ms.estimateRanges(); end
                margins(i) = mlb;
                if mlb <= 0, rob = 2; end
            end
        end

        function ce = falsifyBox(M, lb, ub, label, nTry)
            % random search for a misclassification inside the NORMALISED box.
            ce = false; lab = label + 1;
            for t = 1:nTry
                xn = lb + (ub - lb).*rand(numel(lb),1);   % normalised sample
                xt = ViTReach.unflatten(xn);
                logit = ViTReach.evaluate(M, xt);
                [~, pred] = max(logit);
                if pred ~= lab, ce = true; return; end
            end
        end

        function xt = unflatten(xn)
            % inverse of the (c,h,w) flattening used by epsBox/patchEmbed
            xt = zeros(3,32,32);
            for c = 1:3
                for h = 1:32
                    for w = 1:32
                        xt(c,h,w) = xn((c-1)*1024 + (h-1)*32 + (w-1) + 1);
                    end
                end
            end
        end

        % ---- reach helpers ----
        function Y = bnStar(X, scale, shift, N)
            d = kron(scale(:), ones(N,1));
            b = kron(shift(:), ones(N,1));
            Y = X.affineMap(diag(d), b);
        end

        function R = reluStar(S, method, forced)
            % Sound ReLU on a Star (optionally with forced phases for BaB).
            %   'fast' (default): symbolic triangle relaxation using cheap ESTIMATE
            %     neuron bounds -- no LP, preserves the predicate prefix. This is
            %     what makes the forward reach run in seconds rather than minutes.
            %   otherwise: delegate to NNV's PosLin func (LP-tight but slow).
            %   forced: K x 2 [neuron, phase] -- those neurons are made EXACT
            %     (active: y=x with x>=0 added; inactive: y=0 with x<=0 added). Sound
            %     for the corresponding input-set branch; tighter (no triangle slack).
            if nargin < 3, forced = zeros(0,2); end
            if nargin < 2 || isempty(method) || strcmpi(method, 'fast')
                R = ViTReach.reluStarApproxFast(S, forced);
            else
                R = PosLin.reach(S, method);  R = R(1);   % (forced phases ignored on the LP path)
            end
        end

        function R = reluStarApproxFast(S, forced)
            % approx-star ReLU with estimate (box) neuron bounds; sound triangle,
            % plus exact handling of any 'forced' neurons (BaB phase splits).
            if nargin < 2, forced = zeros(0,2); end
            [l, u] = S.estimateRanges();
            n = S.dim; nVar = S.nVar;
            c = S.V(:, 1); G = S.V(:, 2:end);

            isForced = false(n,1); phaseOf = zeros(n,1);
            if ~isempty(forced)
                isForced(forced(:,1)) = true; phaseOf(forced(:,1)) = forced(:,2);
            end
            unstable = find(l < 0 & u > 0 & ~isForced);
            posn = find(l >= 0 & ~isForced);
            actF = find(isForced & phaseOf > 0);              % forced active -> identity
            % forced inactive and negative-stable -> 0 (rows left zero)
            nr = numel(unstable);

            newV = zeros(n, nVar + 1 + nr);
            idnt = [posn; actF];
            newV(idnt, 1) = c(idnt);
            if nVar > 0, newV(idnt, 2:nVar+1) = G(idnt, :); end
            for k = 1:nr
                newV(unstable(k), nVar+1+k) = 1;              % unstable: y = fresh pred
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
                xc = c(i); xg = G(i, :); col = nVar + k;
                lam = ui / (ui - li);
                rows(3*k-2, col) = -1;                               rd(3*k-2) = 0;          % a>=0
                rows(3*k-1, 1:nVar) = xg; rows(3*k-1, col) = -1;     rd(3*k-1) = -xc;        % a>=x
                rows(3*k,   1:nVar) = -lam*xg; rows(3*k, col) = 1;   rd(3*k) = lam*(xc - li);% a<=line
                pub(col) = ui;
            end
            % forced half-space constraints (restrict alpha to the branch)
            fN = find(isForced);
            frows = zeros(numel(fN), nVar + nr); frd = zeros(numel(fN), 1);
            for k = 1:numel(fN)
                i = fN(k); xc = c(i); xg = G(i, :);
                if phaseOf(i) > 0          % active branch: x_i >= 0  -> -xg*alpha <= xc
                    frows(k, 1:nVar) = -xg; frd(k) = xc;
                else                       % inactive branch: x_i <= 0 ->  xg*alpha <= -xc
                    frows(k, 1:nVar) =  xg; frd(k) = -xc;
                end
            end
            R = Star(newV, [Cold; rows; frows], [dold; rd; frd], plb, pub);
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
