classdef ViTCrown
    % ViTCrown - NNV-native optimized linear-bound reachability for the VNN-COMP 2023
    % ViT. This is NNV's OWN method: the same optimizations NNV's gpu_bab already
    % implements for ReLU nets (gpu_bab_crown_* = backward linear bounds; alpha-
    % relaxation PGA; beta/branch-and-bound), here EXTENDED to the attention-specific
    % layers (bilinear QK^T / A*V and softmax). It is NOT a port of auto_LiRPA's
    % "alpha-CROWN" product - no external dependency; auto_LiRPA is a math reference
    % only. The key tightening lever for the ViT is the alpha-relaxation on the
    % ATTENTION layers (the McCormick bilinear planes), not the FF ReLU, plus score
    % branch-and-bound. Soundness is the only hard constraint (worst case: unknown).
    % Plan: gpu_bab/VIT_ALPHABETA_PLAN.md.
    %
    % The ViT lowers to a flat op DAG that is EXACT affine everywhere except: the FF
    % ReLU, the two bilinear matmuls (QK^T, A*V), and the softmax. This file is the
    % lowering (toOps) + the concrete evaluator (evalOps, for the orientation gate)
    % first; the forward IBP and the backward CROWN bound follow.
    %
    % Op struct fields: .type, .in (source op indices, 1-based into the ops cell),
    % and type-specific data. The value of every op is a column vector; matrix-valued
    % ops carry .mat=[r c] (column-major flatten, reshape(v,[r c]) recovers it).
    %   types: 'input' | 'affine'(.W,.b) | 'relu' | 'add' |
    %          'bmatmul'(.mode 'abt'|'ab', .ra,.ca,.rb,.cb) | 'softmax'(.mat=[N n])
    %
    % Author: NNV Team (ViT track). Date 2026.

    methods (Static)

        % ================= lowering =================
        function ops = toOps(M)
            % Lower the loaded ViTReach model M to a CROWN op DAG (e-major token
            % convention, mirroring ViTReach.reach exactly so the orientation gate
            % can validate it against ViTReach.evaluate).
            E = M.dim; H = M.heads; D = M.dim_head; N = M.N;
            ops = {};
            ops{end+1} = struct('type','input','in',[],'dim',3*32*32,'mat',[]);
            inIdx = 1;

            [Wpe, bpe] = ViTReach.patchEmbed(M);                 % token state dim N*E
            ops{end+1} = ViTCrown.affineOp(inIdx, Wpe, bpe, [N E]);
            x = numel(ops);                                      % current residual-stream op

            for i = 1:M.depth
                blk = M.blocks(i);
                % ---- attention sublayer ----
                xn = ViTCrown.bnOp(ops, x, blk.bn1_scale, blk.bn1_shift, N, E);
                ops{end+1} = xn; xnIdx = numel(ops);
                % Q,K,V per-token linears (e-major: kron(W', I_N))
                Q = ViTCrown.perTokenOp(ops, xnIdx, blk.Mq, blk.bq, N, E); ops{end+1}=Q; qIdx=numel(ops);
                K = ViTCrown.perTokenOp(ops, xnIdx, blk.Mk, blk.bk, N, E); ops{end+1}=K; kIdx=numel(ops);
                V = ViTCrown.perTokenOp(ops, xnIdx, blk.Mv, blk.bv, N, E); ops{end+1}=V; vIdx=numel(ops);
                headOuts = zeros(1,H);
                for h = 1:H
                    cols = (h-1)*D + (1:D);
                    % slice head h's channels (e-major: state idx (e-1)*N+n)
                    Sh = ViTCrown.sliceOp(qIdx, h, D, N, E); ops{end+1}=Sh; qhIdx=numel(ops);
                    Kh = ViTCrown.sliceOp(kIdx, h, D, N, E); ops{end+1}=Kh; khIdx=numel(ops);
                    Vh = ViTCrown.sliceOp(vIdx, h, D, N, E); ops{end+1}=Vh; vhIdx=numel(ops);
                    % scores = scale * Qh * Kh'  (bmatmul mode 'abt'), then scale via affine
                    ops{end+1} = struct('type','bmatmul','in',[qhIdx khIdx],'mode','abt', ...
                        'ra',N,'ca',D,'rb',N,'cb',D,'dim',N*N,'mat',[N N]); scIdx=numel(ops);
                    ops{end+1} = ViTCrown.affineOp(scIdx, M.scale*speye(N*N), zeros(N*N,1), [N N]); scsIdx=numel(ops);
                    % softmax over keys, DECOMPOSED so the backward bound carries a
                    % coefficient to the score (exp -> rowsum -> reciprocal -> broadcast
                    % -> product). This is what lets the QK correlation flow backward
                    % (vs the old constant-box softmax that zeroed the score coeff).
                    ops{end+1} = struct('type','exp','in',scsIdx,'dim',N*N,'mat',[N N]); expIdx=numel(ops);
                    Wsum = sparse(N, N*N);                       % row i sums keys k=1..N of query i
                    for ii=1:N, Wsum(ii, ((0:N-1)*N)+ii) = 1; end
                    ops{end+1} = ViTCrown.affineOp(expIdx, Wsum, zeros(N,1), [N 1]); sumIdx=numel(ops);
                    ops{end+1} = struct('type','reciprocal','in',sumIdx,'dim',N,'mat',[N 1]); recIdx=numel(ops);
                    Wexp = sparse(N*N, N);                       % broadcast R[i] to all keys
                    for kk=1:N, for ii=1:N, Wexp((kk-1)*N+ii, ii)=1; end, end
                    ops{end+1} = ViTCrown.affineOp(recIdx, Wexp, zeros(N*N,1), [N N]); rexpIdx=numel(ops);
                    ops{end+1} = struct('type','eprod','in',[expIdx rexpIdx],'dim',N*N,'mat',[N N]); aIdx=numel(ops);
                    % O_h = A * Vh  (bmatmul mode 'ab')
                    ops{end+1} = struct('type','bmatmul','in',[aIdx vhIdx],'mode','ab', ...
                        'ra',N,'ca',N,'rb',N,'cb',D,'dim',N*D,'mat',[N D]); ohIdx=numel(ops);
                    headOuts(h) = ohIdx;
                end
                % concat heads (e-major) -> dim N*E, then out-projection
                ops{end+1} = ViTCrown.headConcatOp(headOuts, N, D, H, E); oIdx=numel(ops);
                AO = ViTCrown.perTokenOp(ops, oIdx, blk.Mo, blk.bo, N, E); ops{end+1}=AO; aoIdx=numel(ops);
                % residual x = x + AO
                ops{end+1} = struct('type','add','in',[x aoIdx],'dim',N*E,'mat',[N E]); x=numel(ops);

                % ---- feed-forward sublayer ----
                xn2 = ViTCrown.bnOp(ops, x, blk.bn2_scale, blk.bn2_shift, N, E);
                ops{end+1}=xn2; xn2Idx=numel(ops);
                F1 = ViTCrown.perTokenOp(ops, xn2Idx, blk.ff1_W, blk.ff1_b, N, E); ops{end+1}=F1; f1Idx=numel(ops);
                ops{end+1} = struct('type','relu','in',f1Idx,'dim',ops{f1Idx}.dim,'mat',ops{f1Idx}.mat); rIdx=numel(ops);
                F2 = ViTCrown.perTokenOp(ops, rIdx, blk.ff2_W, blk.ff2_b, N, E); ops{end+1}=F2; f2Idx=numel(ops);
                ops{end+1} = struct('type','add','in',[x f2Idx],'dim',N*E,'mat',[N E]); x=numel(ops);
            end
            % ReduceMean over tokens -> dim E
            Wm = ViTReach.meanMap(N, E);
            ops{end+1} = ViTCrown.affineOp(x, Wm, zeros(E,1), [E 1]); zIdx=numel(ops);
            % head BN (per-channel affine) + head linear
            ops{end+1} = ViTCrown.affineOp(zIdx, diag(M.head_bn_scale), M.head_bn_shift, [E 1]); zbIdx=numel(ops);
            ops{end+1} = ViTCrown.affineOp(zbIdx, M.head_W, M.head_b, [10 1]);   % logits
        end

        % ================= concrete evaluator (orientation gate) =================
        function [logits, val] = evalOps(ops, xvec)
            val = cell(numel(ops),1);
            for k = 1:numel(ops)
                op = ops{k};
                switch op.type
                    case 'input',   val{k} = xvec(:);
                    case 'affine',  val{k} = op.W * val{op.in} + op.b;
                    case 'relu',    val{k} = max(val{op.in}, 0);
                    case 'add',     val{k} = val{op.in(1)} + val{op.in(2)};
                    case 'bmatmul'
                        A = reshape(val{op.in(1)}, [op.ra op.ca]);
                        B = reshape(val{op.in(2)}, [op.rb op.cb]);
                        if strcmp(op.mode,'abt'), Cm = A*B'; else, Cm = A*B; end
                        val{k} = Cm(:);
                    case 'softmax'
                        S = reshape(val{op.in}, op.mat);
                        m = max(S,[],2); ex = exp(S-m); A = ex./sum(ex,2);
                        val{k} = A(:);
                    case 'exp',        val{k} = exp(val{op.in});
                    case 'reciprocal', val{k} = 1 ./ val{op.in};
                    case 'eprod',      val{k} = val{op.in(1)} .* val{op.in(2)};
                    case 'concat'
                        % e-major: head h occupies state block (h-1)*N*D+(1:N*D),
                        % so concat is a straight vertical stack of head outputs.
                        val{k} = cell2mat(arrayfun(@(j) val{j}, op.in(:), 'uni', 0));
                    otherwise, error('ViTCrown:evalOps','unknown op %s', op.type);
                end
            end
            logits = val{end};
        end

        % ================= instance verification =================
        function [robust, margins, info] = verifyInstance(ops, M, img01, label, opt)
            % Sound CROWN verification of one benchmark instance (argmax preserved at
            % L-inf eps=1/255). img01 is raw [3,32,32] in [0,1]; label is 0-indexed.
            % opt.eps (default 1/255), opt.alpha (run alpha-opt), opt.alphaIter.
            if nargin < 5, opt = struct(); end
            eps = ViTCrown.getf(opt,'eps',1/255);
            useAlpha = ViTCrown.getf(opt,'alpha',false);
            [lb, ub] = ViTReach.epsBox(M, img01, eps);
            t = label + 1; others = setdiff(1:10, t);
            C = zeros(9,10); for r=1:9, C(r,t)=1; C(r,others(r))=-1; end
            [cl, cu] = ViTCrown.forwardIBP(ops, lb, ub);
            mmaps = ViTCrown.precomputeMaps(ops, cl, cu);
            if useAlpha
                ao = struct('nIter',ViTCrown.getf(opt,'alphaIter',40),'lr',ViTCrown.getf(opt,'lr',0.1));
                margins = ViTCrown.optimizeAlpha(ops, lb, ub, cl, cu, C, ao);
            else
                margins = ViTCrown.backwardCROWN(ops, lb, ub, cl, cu, C, [], mmaps);
            end
            robust = all(margins > 0);
            info = struct('lb',lb,'ub',ub,'C',C);
        end

        % ================= forward IBP (M1) =================
        function [cl, cu] = forwardIBP(ops, lb, ub, ov)
            % Sound per-op interval bounds over the input box [lb,ub]. cl{k},cu{k} are
            % column vectors in op k's value layout. bmatmul = Rump intervalMatMul,
            % softmax = exact correlatedRowSoftmaxBounds (the per-target IBP optimum).
            % ov (optional) = cell(n,1) of [lo hi] range overrides (GenBaB score split):
            % after computing op k its box is intersected with ov{k} (SOUND: case-split
            % on the score value covers the original range).
            n = numel(ops); cl = cell(n,1); cu = cell(n,1);
            if nargin < 4, ov = []; end
            for k = 1:n
                op = ops{k};
                switch op.type
                    case 'input'
                        cl{k} = lb(:); cu{k} = ub(:);
                    case 'affine'
                        Wp = max(op.W,0); Wn = min(op.W,0);
                        cl{k} = Wp*cl{op.in} + Wn*cu{op.in} + op.b;
                        cu{k} = Wp*cu{op.in} + Wn*cl{op.in} + op.b;
                    case 'relu'
                        cl{k} = max(cl{op.in},0); cu{k} = max(cu{op.in},0);
                    case 'add'
                        cl{k} = cl{op.in(1)} + cl{op.in(2)};
                        cu{k} = cu{op.in(1)} + cu{op.in(2)};
                    case 'bmatmul'
                        Al = reshape(cl{op.in(1)},[op.ra op.ca]); Au = reshape(cu{op.in(1)},[op.ra op.ca]);
                        Bl = reshape(cl{op.in(2)},[op.rb op.cb]); Bu = reshape(cu{op.in(2)},[op.rb op.cb]);
                        if strcmp(op.mode,'abt')   % A*B'  (QK^T)
                            [Cl,Cu] = SoftmaxAttn.intervalMatMul(Al,Au,Bl',Bu');
                        else                        % A*B   (A*V)
                            [Cl,Cu] = SoftmaxAttn.intervalMatMul(Al,Au,Bl,Bu);
                        end
                        cl{k} = Cl(:); cu{k} = Cu(:);
                    case 'softmax'
                        Slo = reshape(cl{op.in},op.mat); Shi = reshape(cu{op.in},op.mat);
                        [alb,aub] = SoftmaxAttn.correlatedRowSoftmaxBounds(Slo,Shi);
                        cl{k} = alb(:); cu{k} = aub(:);
                    case 'exp'
                        cl{k} = exp(cl{op.in}); cu{k} = exp(cu{op.in});
                    case 'reciprocal'
                        lo = max(cl{op.in}, 1e-12);             % T = sum exp > 0 (FP floor)
                        cl{k} = 1./cu{op.in}; cu{k} = 1./lo;
                    case 'eprod'
                        El=cl{op.in(1)}; Eu=cu{op.in(1)}; Rl=cl{op.in(2)}; Ru=cu{op.in(2)};
                        c1=El.*Rl; c2=El.*Ru; c3=Eu.*Rl; c4=Eu.*Ru;
                        cl{k} = min(min(c1,c2),min(c3,c4)); cu{k} = max(max(c1,c2),max(c3,c4));
                    case 'concat'
                        cl{k} = cell2mat(arrayfun(@(j) cl{j}, op.in(:),'uni',0));
                        cu{k} = cell2mat(arrayfun(@(j) cu{j}, op.in(:),'uni',0));
                    otherwise, error('ViTCrown:forwardIBP','unknown op %s',op.type);
                end
                if ~isempty(ov) && ~isempty(ov{k})         % GenBaB score range split
                    cl{k} = max(cl{k}, ov{k}(:,1));
                    cu{k} = min(cu{k}, ov{k}(:,2));
                end
            end
        end

        % ================= backward CROWN lower bound (M1b) =================
        function [lbnd, Lin, bL] = backwardCROWN(ops, lb, ub, cl, cu, C, alpha, mmaps, battn, startOp)
            % Lower bound of C*logits over the input box [lb,ub], via backward CROWN.
            % Sound relaxations: affine/add/concat exact; relu triangle (lower slope
            % alpha, upper chord); bmatmul = per-term McCormick (D-reduced) adjoint;
            % softmax = constant IBP box (zero coeff to score). alpha (optional) is a
            % cell: alpha{k} = lower-slope vector for relu op k (default area-adaptive).
            % mmaps (optional) = precomputed bmatmul McCormick maps (precomputeMaps).
            n = numel(ops); S = size(C,1);
            if nargin < 10 || isempty(startOp), startOp = n; end
            Lam = cell(n,1);                 % Lam{k}: [S x dim(op k)] lower-bound coeff
            Lam{startOp} = C; bL = zeros(S,1);   % spec on startOp's output (n = logits)
            if nargin < 7 || isempty(alpha), alpha = cell(n,1); end
            if nargin < 8 || isempty(mmaps), mmaps = ViTCrown.precomputeMaps(ops, cl, cu); end
            if nargin < 9 || isempty(battn), battn = cell(n,1); end
            for k = n:-1:1
                L = Lam{k};
                if isempty(L), continue; end
                op = ops{k};
                switch op.type
                    case 'input'
                        % leaf: keep L as Lin
                    case 'affine'
                        Lam{op.in} = ViTCrown.addco(Lam{op.in}, L*op.W);
                        bL = bL + L*op.b;
                    case 'add'
                        Lam{op.in(1)} = ViTCrown.addco(Lam{op.in(1)}, L);
                        Lam{op.in(2)} = ViTCrown.addco(Lam{op.in(2)}, L);
                    case 'concat'
                        off = 0;
                        for j = op.in(:)'
                            dj = numel(cl{j});
                            Lam{j} = ViTCrown.addco(Lam{j}, L(:, off+(1:dj)));
                            off = off + dj;
                        end
                    case 'relu'
                        l = cl{op.in}; u = cu{op.in};
                        sp = double(l >= 0); sn = double(u <= 0);     % stable +/-
                        un = 1 - sp - sn;                             % unstable mask
                        den = (u - l); den(den==0) = 1;
                        slope_u = u./den;
                        su = sp + un.*slope_u;                        % upper-chord slope
                        ci = un.*(-slope_u.*l);                       % upper-chord intercept
                        if isempty(alpha{k})
                            a = double(u > -l);                       % area-adaptive default
                        else
                            a = alpha{k}(:);                          % may be dlarray (M5)
                        end
                        sL = sp + un.*a;                              % lower slope (masked, dlarray-safe)
                        Lp = max(L,0); Ln = min(L,0);
                        Lam{op.in} = ViTCrown.addco(Lam{op.in}, Lp.*sL.' + Ln.*su.');
                        bL = bL + Ln*ci;
                    case 'softmax'
                        % constant box: y in [cl,cu]; zero coeff to score
                        Lp = max(L,0); Ln = min(L,0);
                        bL = bL + Lp*cl{k} + Ln*cu{k};
                    case 'exp'
                        % convex: lower tangent at contact m, upper chord. Carries a
                        % real coefficient to the score (the QK correlation lever).
                        l = cl{op.in}; u = cu{op.in}; m = (l+u)/2;
                        slo = exp(m); ilo = exp(m).*(1-m);              % lower: slo*x+ilo
                        du = u-l; du(du==0)=1; sup = (exp(u)-exp(l))./du; iup = exp(l)-sup.*l;
                        Lp = max(L,0); Ln = min(L,0);
                        Lam{op.in} = ViTCrown.addco(Lam{op.in}, Lp.*slo.' + Ln.*sup.');
                        bL = bL + Lp*ilo + Ln*iup;
                    case 'reciprocal'
                        % 1/x convex on x>0: lower tangent at m, upper chord
                        l = max(cl{op.in},1e-12); u = max(cu{op.in},1e-12); m = (l+u)/2;
                        slo = -1./(m.^2); ilo = 2./m;                   % lower tangent
                        sup = -1./(l.*u); iup = 1./l - sup.*l;          % upper chord
                        Lp = max(L,0); Ln = min(L,0);
                        Lam{op.in} = ViTCrown.addco(Lam{op.in}, Lp.*slo.' + Ln.*sup.');
                        bL = bL + Lp*ilo + Ln*iup;
                    case 'eprod'
                        % elementwise product E.*R, per-element McCormick (default pick)
                        El=cl{op.in(1)}; Eu=cu{op.in(1)}; Rl=cl{op.in(2)}; Ru=cu{op.in(2)};
                        [PL,PU,mid] = ViTCrown.mcCorners(El,Eu,Rl,Ru);
                        aL = double(mid.L1>=mid.L2); aU = double(mid.U1<=mid.U2);
                        apL=aL.*PL.ap1+(1-aL).*PL.ap2; aqL=aL.*PL.aq1+(1-aL).*PL.aq2; a0L=aL.*PL.a01+(1-aL).*PL.a02;
                        apU=aU.*PU.ap1+(1-aU).*PU.ap2; aqU=aU.*PU.aq1+(1-aU).*PU.aq2; a0U=aU.*PU.a01+(1-aU).*PU.a02;
                        Lp = max(L,0); Ln = min(L,0);
                        Lam{op.in(1)} = ViTCrown.addco(Lam{op.in(1)}, Lp.*apL.' + Ln.*apU.');
                        Lam{op.in(2)} = ViTCrown.addco(Lam{op.in(2)}, Lp.*aqL.' + Ln.*aqU.');
                        bL = bL + Lp*a0L + Ln*a0U;
                    case 'bmatmul'
                        mp = mmaps{k};
                        if ~isempty(battn{k}), aL = battn{k}.aL; aU = battn{k}.aU; else, aL=[]; aU=[]; end
                        [PAL,PAU,PBL,PBU,p0L,p0U] = ViTCrown.bmmCombine(mp, aL, aU);
                        Lp = max(L,0); Ln = min(L,0);
                        Lam{op.in(1)} = ViTCrown.addco(Lam{op.in(1)}, Lp*PAL + Ln*PAU);
                        Lam{op.in(2)} = ViTCrown.addco(Lam{op.in(2)}, Lp*PBL + Ln*PBU);
                        bL = bL + Lp*p0L + Ln*p0U;
                    otherwise, error('ViTCrown:backwardCROWN','unknown op %s',op.type);
                end
            end
            Lin = Lam{1};
            lbnd = max(Lin,0)*lb + min(Lin,0)*ub + bL;   % min over box
        end

        % ============= CROWN intermediate bounds (alpha,beta-CROWN style) ====
        function [lo, hi] = crownBounds(ops, lb, ub, cl, cu, targetOp, mmaps)
            % Tighter [lo,hi] for op targetOp's output via a backward CROWN pass to
            % the input (vs the looser forward-IBP). lo = min, hi = -min(-x).
            d = ops{targetOp}.dim; I = eye(d);
            lo = ViTCrown.backwardCROWN(ops, lb, ub, cl, cu,  I, [], mmaps, [], targetOp);
            hi = -ViTCrown.backwardCROWN(ops, lb, ub, cl, cu, -I, [], mmaps, [], targetOp);
        end

        function [cl, cu, mmaps] = refineScores(ops, lb, ub)
            % Replace the IBP score boxes (exp inputs) with tighter CROWN bounds, then
            % re-propagate IBP downstream (so the exp/A/A*V relaxations see the tighter
            % scores). This is the key alpha,beta-CROWN tightness lever for attention.
            [cl,cu] = ViTCrown.forwardIBP(ops, lb, ub);
            mmaps = ViTCrown.precomputeMaps(ops, cl, cu);
            expIdx = find(cellfun(@(o) strcmp(o.type,'exp'), ops));
            ov = cell(numel(ops),1);
            for e = expIdx(:)'
                sc = ops{e}.in;
                [lo,hi] = ViTCrown.crownBounds(ops, lb, ub, cl, cu, sc, mmaps);
                ov{sc} = [max(cl{sc},lo(:)), min(cu{sc},hi(:))];
            end
            [cl,cu] = ViTCrown.forwardIBP(ops, lb, ub, ov);
            mmaps = ViTCrown.precomputeMaps(ops, cl, cu);
        end

        function [cl, cu, mmaps] = refineBounds(ops, lb, ub, iters, maxDim)
            % Replace IBP boxes with CROWN intermediate bounds for ALL nonlinearity
            % inputs (relu/exp/reciprocal/eprod/bmatmul operands), sequentially in
            % topological order so later refinements use earlier ones, re-propagating
            % IBP each pass. This is the alpha,beta-CROWN tightness driver - with it the
            % backward bound approaches the LP optimum. maxDim caps which ops to refine
            % (skip very wide ones for speed).
            if nargin < 4 || isempty(iters), iters = 1; end
            if nargin < 5 || isempty(maxDim), maxDim = inf; end
            [cl,cu] = ViTCrown.forwardIBP(ops, lb, ub);
            mmaps = ViTCrown.precomputeMaps(ops, cl, cu);
            nl = {'relu','exp','reciprocal','eprod','bmatmul'};
            targets = [];
            for k = 1:numel(ops)
                if any(strcmp(ops{k}.type, nl)), targets = [targets, ops{k}.in(:)']; end %#ok<AGROW>
            end
            targets = unique(targets);
            targets = targets(arrayfun(@(t) ops{t}.dim <= maxDim, targets));
            for it = 1:iters
                ov = cell(numel(ops),1);
                for tt = sort(targets)
                    [lo,hi] = ViTCrown.crownBounds(ops, lb, ub, cl, cu, tt, mmaps);
                    ov{tt} = [max(cl{tt},lo(:)), min(cu{tt},hi(:))];
                end
                [cl,cu] = ViTCrown.forwardIBP(ops, lb, ub, ov);
                mmaps = ViTCrown.precomputeMaps(ops, cl, cu);
            end
        end

        function mmaps = precomputeMaps(ops, cl, cu)
            % Precompute, for every bmatmul op, the FIXED sparse McCormick adjoint
            % maps (plane selection depends only on the box, not on the spec or alpha).
            % Then the backward pass is pure matmuls: GA = Lp*PAL + Ln*PAU, etc.
            mmaps = cell(numel(ops),1);
            for k = 1:numel(ops)
                if strcmp(ops{k}.type,'bmatmul')
                    mmaps{k} = ViTCrown.bmmMaps(ops{k}, cl, cu);
                end
            end
        end

        function mp = bmmMaps(op, cl, cu)
            % Build the TWO McCormick corner maps for each of the lower/upper envelopes
            % of a reducing bilinear matmul. A sound linear under-estimator is any
            % convex combo aL*L1+(1-aL)*L2 of the two McCormick lower planes (likewise
            % aU for upper) -> aL,aU in [0,1] are the ATTENTION-LAYER alpha (optimized
            % in optimizeAlpha). Stored dense (dimY<=289) so they flow through dlarray.
            clA = reshape(cl{op.in(1)},[op.ra op.ca]); cuA = reshape(cu{op.in(1)},[op.ra op.ca]);
            clB = reshape(cl{op.in(2)},[op.rb op.cb]); cuB = reshape(cu{op.in(2)},[op.rb op.cb]);
            dimY = op.dim; dimA = op.ra*op.ca; dimB = op.rb*op.cb; N = op.ra;
            rA=[];cA=[]; AL1=[];AL2=[];AU1=[];AU2=[];
            rB=[];cB=[]; BL1=[];BL2=[];BU1=[];BU2=[];
            c0L1=zeros(dimY,1); c0L2=zeros(dimY,1); c0U1=zeros(dimY,1); c0U2=zeros(dimY,1);
            aL0=zeros(dimY,1); aU0=zeros(dimY,1);
            if strcmp(op.mode,'abt'), Iout=op.ra; Jout=op.rb; else, Iout=op.cb; Jout=op.ra; end
            for a1 = 1:Iout
                for a2 = 1:Jout
                    if strcmp(op.mode,'abt')   % Y(i,j)=sum_d A(i,d)B(j,d)
                        i=a1; j=a2; m=(j-1)*N+i;
                        pl=clA(i,:).'; pu=cuA(i,:).'; idxA=((0:op.ca-1)*N+i).';
                        ql=clB(j,:).'; qu=cuB(j,:).'; idxB=((0:op.cb-1)*N+j).';
                    else                        % Y(i,e)=sum_n A(i,n)B(n,e)
                        e=a1; i=a2; m=(e-1)*N+i;
                        pl=clA(i,:).'; pu=cuA(i,:).'; idxA=((0:op.ca-1)*N+i).';
                        ql=clB(:,e);   qu=cuB(:,e);   idxB=((e-1)*N+(1:op.rb)).';
                    end
                    [PL,PU,mid] = ViTCrown.mcCorners(pl,pu,ql,qu);
                    nA=numel(idxA); nB=numel(idxB);
                    rA=[rA;m*ones(nA,1)]; cA=[cA;idxA]; AL1=[AL1;PL.ap1]; AL2=[AL2;PL.ap2]; AU1=[AU1;PU.ap1]; AU2=[AU2;PU.ap2];
                    rB=[rB;m*ones(nB,1)]; cB=[cB;idxB]; BL1=[BL1;PL.aq1]; BL2=[BL2;PL.aq2]; BU1=[BU1;PU.aq1]; BU2=[BU2;PU.aq2];
                    c0L1(m)=sum(PL.a01); c0L2(m)=sum(PL.a02); c0U1(m)=sum(PU.a01); c0U2(m)=sum(PU.a02);
                    aL0(m)=double(sum(mid.L1)>=sum(mid.L2));
                    aU0(m)=double(sum(mid.U1)<=sum(mid.U2));
                end
            end
            mp.AL1=full(sparse(rA,cA,AL1,dimY,dimA)); mp.AL2=full(sparse(rA,cA,AL2,dimY,dimA));
            mp.AU1=full(sparse(rA,cA,AU1,dimY,dimA)); mp.AU2=full(sparse(rA,cA,AU2,dimY,dimA));
            mp.BL1=full(sparse(rB,cB,BL1,dimY,dimB)); mp.BL2=full(sparse(rB,cB,BL2,dimY,dimB));
            mp.BU1=full(sparse(rB,cB,BU1,dimY,dimB)); mp.BU2=full(sparse(rB,cB,BU2,dimY,dimB));
            mp.c0L1=c0L1; mp.c0L2=c0L2; mp.c0U1=c0U1; mp.c0U2=c0U2;
            mp.aL0=aL0; mp.aU0=aU0; mp.dimY=dimY;
        end

        function [PAL,PAU,PBL,PBU,p0L,p0U] = bmmCombine(mp, aL, aU)
            % Effective McCormick maps for per-output-element attention-alpha aL,aU in
            % [0,1] (convex combo of two valid planes -> sound for any a).
            if isempty(aL), aL = mp.aL0; end
            if isempty(aU), aU = mp.aU0; end
            PAL = aL.*mp.AL1 + (1-aL).*mp.AL2;  PBL = aL.*mp.BL1 + (1-aL).*mp.BL2;
            PAU = aU.*mp.AU1 + (1-aU).*mp.AU2;  PBU = aU.*mp.BU1 + (1-aU).*mp.BU2;
            p0L = aL.*mp.c0L1 + (1-aL).*mp.c0L2;  p0U = aU.*mp.c0U1 + (1-aU).*mp.c0U2;
        end

        function A = addco(A, B)
            % accumulate coefficient block B into A (A may be [] = zero)
            if isempty(A), A = B; else, A = A + B; end
        end

        function v = getf(s, f, d)
            if isfield(s,f), v = s.(f); else, v = d; end
        end

        % ================= alpha-relaxation optimization =================
        function [lbnd, sol, hist] = optimizeAlpha(ops, lb, ub, cl, cu, C, opt)
            % Projected-gradient ascent (Adam, dlgradient) on the relaxation alpha to
            % maximize the worst margin. Optimizes BOTH the FF-relu lower slope AND -
            % the main lever for the ViT - the ATTENTION-LAYER alpha: the per-output
            % McCormick plane interpolation aL,aU on every bilinear QK^T / A*V op.
            % The optimized bound is RE-EVALUATED on the sound double path, so
            % soundness never depends on the autodiff arithmetic (any alpha in [0,1]
            % is sound). Returns sol = {aReluCell, battnCell} for reuse by BaB.
            if nargin < 7, opt = struct(); end
            nIter = ViTCrown.getf(opt,'nIter',60); lr = ViTCrown.getf(opt,'lr',0.05);
            verbose = ViTCrown.getf(opt,'verbose',false);
            doRelu = ViTCrown.getf(opt,'relu',true); doAttn = ViTCrown.getf(opt,'attn',true);
            mmaps = ViTCrown.precomputeMaps(ops, cl, cu);
            reluIdx = find(cellfun(@(o) strcmp(o.type,'relu'), ops));
            bmIdx   = find(cellfun(@(o) strcmp(o.type,'bmatmul'), ops));
            % variable list with metadata (kind, op index)
            vars = {}; meta = {};
            if doRelu
                for r = 1:numel(reluIdx)
                    inb = ops{reluIdx(r)}.in;
                    vars{end+1} = dlarray(double(cu{inb} > -cl{inb})); meta{end+1} = {'relu', reluIdx(r)}; %#ok<AGROW>
                end
            end
            if doAttn
                for b = 1:numel(bmIdx)
                    mp = mmaps{bmIdx(b)};
                    vars{end+1} = dlarray(mp.aL0); meta{end+1} = {'aL', bmIdx(b)}; %#ok<AGROW>
                    vars{end+1} = dlarray(mp.aU0); meta{end+1} = {'aU', bmIdx(b)}; %#ok<AGROW>
                end
            end
            mAdam = cellfun(@(a) zeros(size(a)), vars, 'uni', 0); vAdam = mAdam;
            b1=0.9; b2=0.999; eps=1e-8; best=-inf; bestV=vars; hist=zeros(nIter,1);
            for it = 1:nIter
                [obj, grads] = dlfeval(@ViTCrown.alphaObj, vars, meta, ops, lb, ub, cl, cu, C, mmaps);
                ov = double(gather(extractdata(obj))); hist(it) = ov;
                if ov > best, best = ov; bestV = vars; end
                for v = 1:numel(vars)
                    g = grads{v};
                    mAdam{v} = b1*mAdam{v} + (1-b1)*g; vAdam{v} = b2*vAdam{v} + (1-b2)*g.^2;
                    mhat = mAdam{v}/(1-b1^it); vhat = vAdam{v}/(1-b2^it);
                    vars{v} = max(min(vars{v} + lr * mhat./(sqrt(vhat)+eps), 1), 0);   % ascent + project
                end
                if verbose && (it<=3 || mod(it,10)==0)
                    fprintf('   alpha-opt it %3d: min-margin = %.4f\n', it, ov);
                end
            end
            [aReluCell, battn] = ViTCrown.buildAlpha(bestV, meta, numel(ops), true);
            lbnd = ViTCrown.backwardCROWN(ops, lb, ub, cl, cu, C, aReluCell, mmaps, battn);
            sol = struct('aRelu',{aReluCell},'battn',{battn});
        end

        function [aReluCell, battn] = buildAlpha(vars, meta, n, toDouble)
            % assemble backwardCROWN's relu-alpha cell and bmatmul-attn (aL/aU) cell
            aReluCell = cell(n,1); battn = cell(n,1);
            for v = 1:numel(vars)
                m = meta{v}; val = vars{v};
                if toDouble, val = double(gather(extractdata(val))); end
                switch m{1}
                    case 'relu', aReluCell{m{2}} = val;
                    case 'aL', if isempty(battn{m{2}}), battn{m{2}}=struct('aL',[],'aU',[]); end, battn{m{2}}.aL = val;
                    case 'aU', if isempty(battn{m{2}}), battn{m{2}}=struct('aL',[],'aU',[]); end, battn{m{2}}.aU = val;
                end
            end
        end

        function [obj, grads] = alphaObj(vars, meta, ops, lb, ub, cl, cu, C, mmaps)
            [aReluCell, battn] = ViTCrown.buildAlpha(vars, meta, numel(ops), false);
            lbnd = ViTCrown.backwardCROWN(ops, lb, ub, cl, cu, C, aReluCell, mmaps, battn);
            obj = min(lbnd);                       % worst margin (maximize)
            grads = dlgradient(obj, vars);
        end

        function [PL,PU,mid] = mcCorners(pl,pu,ql,qu)
            % The two McCormick under-estimator planes (L1,L2) and over-estimator
            % planes (U1,U2) for t=p*q on [pl,pu]x[ql,qu]. Each as coeffs ap*p+aq*q+a0.
            % mid.* = plane value at the box midpoint (for the default tighter pick).
            pm=(pl+pu)/2; qm=(ql+qu)/2;
            PL.ap1=ql; PL.aq1=pl; PL.a01=-(pl.*ql);          % pl*q+ql*p-pl*ql
            PL.ap2=qu; PL.aq2=pu; PL.a02=-(pu.*qu);          % pu*q+qu*p-pu*qu
            mid.L1 = pl.*qm+ql.*pm-pl.*ql;  mid.L2 = pu.*qm+qu.*pm-pu.*qu;
            PU.ap1=ql; PU.aq1=pu; PU.a01=-(pu.*ql);          % pu*q+ql*p-pu*ql
            PU.ap2=qu; PU.aq2=pl; PU.a02=-(pl.*qu);          % pl*q+qu*p-pl*qu
            mid.U1 = pu.*qm+ql.*pm-pu.*ql;  mid.U2 = pl.*qm+qu.*pm-pl.*qu;
        end

        % ================= op builders =================
        function op = affineOp(inIdx, W, b, mat)
            op = struct('type','affine','in',inIdx,'W',W,'b',b(:),'dim',size(W,1),'mat',mat);
        end
        function op = perTokenOp(ops, inIdx, Mw, b, N, E) %#ok<INUSL>
            % Y = Xt*Mw + b per token; e-major map = kron(Mw', I_N). Mw is [Ein x Eout].
            Eout = size(Mw,2);
            W = kron(Mw', speye(N)); bb = kron(b(:), ones(N,1));
            op = ViTCrown.affineOp(inIdx, W, bb, [N Eout]);
        end
        function op = bnOp(ops, inIdx, scale, shift, N, E) %#ok<INUSL>
            d = kron(scale(:), ones(N,1)); bsh = kron(shift(:), ones(N,1));
            op = ViTCrown.affineOp(inIdx, spdiags(d,0,N*E,N*E), bsh, [N E]);
        end
        function op = sliceOp(inIdx, h, D, N, E) %#ok<INUSD>
            % select head h channels (e-major contiguous block) -> dim N*D
            idx = ((h-1)*D)*N + (1:D*N);
            W = sparse(1:numel(idx), idx, 1, numel(idx), N*E);
            op = ViTCrown.affineOp(inIdx, W, zeros(numel(idx),1), [N D]);
        end
        function op = headConcatOp(headOuts, N, D, H, E)
            % Stack head outputs (each dim N*D) into e-major dim N*E. In e-major, head
            % h occupies the contiguous state block (h-1)*N*D + (1:N*D), so the concat
            % is exactly a vertical stack of the head-output vectors in head order.
            op = struct('type','concat','in',headOuts,'N',N,'D',D,'H',H,'E',E,'dim',N*E,'mat',[N E]);
        end

    end
end
