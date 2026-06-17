function gpu_bab_alpha_dag_parity_test()
% GPU_BAB_ALPHA_DAG_PARITY_TEST  Soundness gate for gpu_bab_crown_alpha_dag.
%   Asserts gpu_bab_crown_alpha_dag(nIter=0) == gpu_bab_crown_spec_dag BOUND-FOR-BOUND on
%   synthetic conv DAGs (conv/relu/avgpool/affine/normaffine/add). At min-area alpha the two
%   backwards are op-for-op identical, so any mismatch is a fused-backward bug -- this MUST pass
%   before alpha_dag is trusted for a 'robust' emit. Also checks (a) nIter>0 runs and is
%   sound-AND-no-worse than min-area, and (b) the rootBounds-reuse path matches too.
%
%   Run: results = runtests('gpu_bab_alpha_dag_parity_test')  OR  gpu_bab_alpha_dag_parity_test
    here = fileparts(mfilename('fullpath'));
    addpath(here);
    rng(7, 'twister');
    prec = 'double';                                  % double: parity should be ~exact
    tol = 1e-9;

    nets = {i_net_chain(), i_net_residual()};
    names = {'chain conv->relu->avgpool->fc', 'normaffine->conv->relu->conv->add->relu->fc'};

    for ni = 1:numel(nets)
        net = nets{ni}; ops = net.ops; reluIdx = net.reluIdx;
        n = net.nIn; nOut = net.nOut; B = 3;
        % random input box + spec
        c  = randn(n, B);
        rad = 0.3 * (0.5 + rand(n, B));
        lb = c - rad; ub = c + rad;
        C  = randn(4, nOut);

        % ---- IBP path: alpha_dag(nIter=0) vs spec_dag ----
        mSpec = gpu_bab_crown_spec_dag(ops, lb, ub, C, prec, {}, []);
        mA0   = gpu_bab_crown_alpha_dag(ops, lb, ub, C, {}, reluIdx, prec, 0, [], []);
        e = max(abs(mSpec(:) - mA0(:)));
        assert(e < tol, '[%s] IBP parity FAILED: max|spec-alpha0| = %.3e (tol %.1e)', names{ni}, e, tol);
        fprintf('[%s] IBP parity OK (max diff %.2e)\n', names{ni}, e);

        % ---- rootBounds-reuse path: build root bounds from the IBP forward, compare ----
        rb = i_root_bounds(ops, lb, ub, reluIdx, prec);
        mSpecR = gpu_bab_crown_spec_dag(ops, lb, ub, C, prec, {}, rb);
        mA0R   = gpu_bab_crown_alpha_dag(ops, lb, ub, C, {}, reluIdx, prec, 0, [], rb);
        eR = max(abs(mSpecR(:) - mA0R(:)));
        assert(eR < tol, '[%s] rootBounds parity FAILED: %.3e', names{ni}, eR);
        fprintf('[%s] rootBounds parity OK (max diff %.2e)\n', names{ni}, eR);

        % ---- nIter>0: runs + sound-AND-no-worse than min-area (per-spec lower bound) ----
        mOpt = gpu_bab_crown_alpha_dag(ops, lb, ub, C, {}, reluIdx, prec, 15, 0.2, []);
        worst0  = min(mA0, [], 1);                    % the certified (worst-spec) margin per node
        worstOpt = min(mOpt, [], 1);
        gain = min(worstOpt - worst0);
        % SOUNDNESS: keep-best guarantees alpha-opt is never worse than min-area.
        assert(all(worstOpt >= worst0 - 1e-7), ...
            '[%s] alpha-opt REGRESSED below min-area (unsound keep-best?): %.3e', names{ni}, gain);
        % ROUTE-B REGRESSION GUARD: with the traceable conv/pool backward, dlgradient must reach
        % alpha and IMPROVE the bound on these nets (gain>0). gain==0 => the autodiff tape broke
        % again (extractdata / indexed-assign crept back in) -> optimizer is a silent no-op.
        assert(gain > 1e-4, ...
            '[%s] alpha-opt did NOT improve (gain %.3e) -- autodiff tape likely severed (Route B).', names{ni}, gain);
        fprintf('[%s] nIter>0 sound + OPTIMIZING OK (worst-margin gain %.3e)\n', names{ni}, gain);
    end

    % ============ BETA soundness (Monte Carlo) — the -150 gate for the BaB path ============
    % With BaB split fixings, alpha_dag optimizes [alpha;beta]. The returned bound MUST be a valid
    % lower bound on min C*f(x) over inputs x in the box that SATISFY the split constraints. Sample
    % such inputs, eval the true margin C*f, and assert bound <= the MC true-min (mB <= true_min <=
    % MC_min). A violation means the beta dual is mis-implemented (wrong sign/neuron) -> false robust.
    for ni = 1:numel(nets)
        net = nets{ni}; ops = net.ops; reluIdx = net.reluIdx; n = net.nIn; nOut = net.nOut;
        c = randn(n,1); rad = 0.4*(0.5+rand(n,1)); lb = c-rad; ub = c+rad; C = randn(3, nOut);
        [fixings, x0] = i_make_fixings(ops, lb, ub, reluIdx, prec);
        mB = gpu_bab_crown_alpha_dag(ops, lb, ub, C, fixings, reluIdx, prec, 20, 0.1, []); % alpha+beta
        mB = min(mB, [], 2);                                  % per-spec bound (B=1)
        nMC = 40000; X = [x0, lb + (ub - lb) .* rand(n, nMC-1)];   % witness x0 -> >=1 feasible sample
        [Y, Z] = i_forward(ops, X, reluIdx);
        ok = i_satisfies(Z, fixings, reluIdx, 0);
        assert(any(ok), '[beta %d] no MC sample satisfied the split fixings (bad test setup)', ni);
        minTrue = min(C * Y(:, ok), [], 2);                  % per-spec true min over feasible samples
        viol = max(mB - minTrue);                            % must be <= 0 (sound)
        assert(viol <= 1e-6, '[beta %d] BETA UNSOUND: bound exceeds MC true-min by %.3e', ni, viol);
        fprintf('[beta soundness %d] OK (bound <= MC true-min; slack %.3e; %d/%d feasible)\n', ...
            ni, min(minTrue - mB), sum(ok), nMC);
    end

    fprintf('\nALL PARITY + BETA-SOUNDNESS TESTS PASSED\n');
end

% =====================================================================================
function rb = i_root_bounds(ops, lb, ub, reluIdx, prec)
% Pack a rootBounds struct (per-relu preL/preU, dim_k x 1) from the IBP forward over the full
% box -- a valid (loose) "tight" bound for the reuse-path parity check.
    nOps = numel(ops);
    cl = cell(nOps+1,1); cu = cell(nOps+1,1);
    cl{1} = cast(lb(:,1), prec); cu{1} = cast(ub(:,1), prec);
    preL = cell(nOps,1); preU = cell(nOps,1);
    for k = 1:nOps
        op = ops{k};
        if strcmp(op.type,'add')
            a = op.inputs(1)+1; b = op.inputs(2)+1;
            cl{k+1} = cl{a}+cl{b}; cu{k+1} = cu{a}+cu{b}; continue;
        end
        s = op.src+1; l = cl{s}; u = cu{s};
        switch op.type
            case 'affine'
                W = op.W; bb = op.b(:); Wp = max(W,0); Wn = min(W,0);
                cl{k+1} = Wp*l + Wn*u + bb; cu{k+1} = Wp*u + Wn*l + bb;
            case 'conv'
                [cl{k+1}, cu{k+1}] = i_conv_ibp_local(op, l, u, prec);
            case 'normaffine'
                sf = i_bcast(op.scale, op.shape); tf = i_bcast(op.shift, op.shape);
                pos = sf>=0;
                cl{k+1} = (sf.*l).*pos + (sf.*u).*(~pos) + tf;
                cu{k+1} = (sf.*u).*pos + (sf.*l).*(~pos) + tf;
            case 'avgpool'
                [cl{k+1}, cu{k+1}] = i_avgpool_ibp_local(op, l, u, prec);
            case 'relu'
                preL{k} = l; preU{k} = u; cl{k+1} = max(l,0); cu{k+1} = max(u,0);
        end
    end
    rb = struct('preL', {preL}, 'preU', {preU});
end

% ---- synthetic nets (hand-built ops; shapes verified consistent with dlconv) -------------
function net = i_net_chain()
    ops = {};
    ops{end+1} = i_conv([4 4 2], 3, 3, [1 1], [1 1 1 1], 0);     % op1: [4 4 2]->[4 4 3]
    ops{end+1} = struct('type','relu','src',1);                  % op2
    ops{end+1} = i_avgpool([4 4 3], [2 2], 2);                   % op3: ->[2 2 3]
    ops{end+1} = i_affine(5, prod([2 2 3]), 3);                  % op4: 12->5
    net = struct('ops', {ops}, 'reluIdx', 2, 'nIn', prod([4 4 2]), 'nOut', 5);
end

function net = i_net_residual()
    ops = {};
    ops{end+1} = i_normaffine([4 4 2], 0);                        % op1
    ops{end+1} = i_conv([4 4 2], 3, 3, [1 1], [1 1 1 1], 1);      % op2: ->[4 4 3]
    ops{end+1} = struct('type','relu','src',2);                  % op3
    ops{end+1} = i_conv([4 4 3], 3, 3, [1 1], [1 1 1 1], 3);      % op4: ->[4 4 3]
    ops{end+1} = struct('type','add','inputs',[3 4],'shape',[4 4 3],'src',3); % op5: op3+op4
    ops{end+1} = struct('type','relu','src',5);                  % op6
    ops{end+1} = i_affine(5, prod([4 4 3]), 6);                  % op7: 48->5
    net = struct('ops', {ops}, 'reluIdx', [3 6], 'nIn', prod([4 4 2]), 'nOut', 5);
end

function op = i_conv(inShape, fh, Cout, stride, pad, src)
    Cin = inShape(3);
    W = 0.3*randn(fh, fh, Cin, Cout); b = 0.1*randn(1,1,Cout);
    Hout = floor((inShape(1)+pad(1)+pad(2)-fh)/stride(1)+1);
    Wout = floor((inShape(2)+pad(3)+pad(4)-fh)/stride(2)+1);
    op = struct('type','conv','W',W,'b',b,'stride',stride,'pad',pad,'dil',[1 1], ...
                'inShape',inShape,'outShape',[Hout Wout Cout],'src',src);
end
function op = i_avgpool(inShape, pool, src)
    Hout = floor(inShape(1)/pool(1)); Wout = floor(inShape(2)/pool(2));
    op = struct('type','avgpool','pool',pool,'stride',pool,'pad',[0 0 0 0], ...
                'inShape',inShape,'outShape',[Hout Wout inShape(3)],'src',src);
end
function op = i_affine(nOut, nIn, src)
    op = struct('type','affine','W',0.3*randn(nOut,nIn),'b',0.1*randn(nOut,1),'src',src);
end
function op = i_normaffine(shape, src)
    C = shape(3);
    op = struct('type','normaffine','scale',0.5+rand(1,1,C),'shift',0.2*randn(1,1,C),'shape',shape,'src',src);
end

% ---- local IBP helpers for i_root_bounds (mirror the engine's) ---------------------------
function [olb, oub] = i_conv_ibp_local(op, lb, ub, prec)
    ish = op.inShape; osh = op.outShape; B = size(lb,2);
    W = cast(op.W, prec); Wp = max(W,0); Wn = min(W,0);
    bb = reshape(cast(op.b(:), prec), [1 1 osh(3)]);
    L4 = dlarray(reshape(lb, [ish(1) ish(2) ish(3) B]), 'SSCB');
    U4 = dlarray(reshape(ub, [ish(1) ish(2) ish(3) B]), 'SSCB');
    pad2 = [op.pad(1) op.pad(3); op.pad(2) op.pad(4)];
    args = {'Stride', op.stride, 'Padding', pad2, 'DilationFactor', op.dil};
    Lo = dlconv(L4, Wp, bb, args{:}) + dlconv(U4, Wn, 0, args{:});
    Hi = dlconv(U4, Wp, bb, args{:}) + dlconv(L4, Wn, 0, args{:});
    olb = reshape(extractdata(Lo), [prod(osh) B]); oub = reshape(extractdata(Hi), [prod(osh) B]);
end
function [olb, oub] = i_avgpool_ibp_local(op, lb, ub, prec)
    ish = op.inShape; osh = op.outShape; B = size(lb,2); kh = op.pool(1); kw = op.pool(2);
    L4 = reshape(cast(lb,prec), [ish(1) ish(2) ish(3) B]); U4 = reshape(cast(ub,prec), [ish(1) ish(2) ish(3) B]);
    olb = zeros([osh prod(1) B]); olb = reshape(i_pm(L4,osh,kh,kw,op.stride), [prod(osh) B]);
    oub = reshape(i_pm(U4,osh,kh,kw,op.stride), [prod(osh) B]);
end
function Y = i_pm(X, osh, kh, kw, stride)
    B = size(X,4); Y = zeros([osh(1) osh(2) osh(3) B], 'like', X);
    for oh=1:osh(1), for ow=1:osh(2)
        rh=(oh-1)*stride(1)+(1:kh); rw=(ow-1)*stride(2)+(1:kw);
        Y(oh,ow,:,:)=mean(mean(X(rh,rw,:,:),1),2);
    end, end
end
function v = i_bcast(x, sh)
    v = reshape(zeros([sh(1) sh(2) sh(3)]) + x, [], 1);
end

% ---- beta-soundness helpers: split fixings + exact forward eval + constraint satisfaction -----
function [fixings, x0] = i_make_fixings(ops, lb, ub, reluIdx, prec)
% Fix up to 4 unstable neurons of the FIRST relu layer to a WITNESS point x0's actual signs, so the
% split region is guaranteed non-empty (x0 satisfies it). IBP-unstable != forward-straddling, so
% fixing by a real eval (not by sign alternation) is what keeps the constraints feasible.
    rb = i_root_bounds(ops, lb, ub, reluIdx, prec);
    nOps = numel(ops); fixings = cell(nOps,1);
    k = reluIdx(1); l = rb.preL{k};
    u = rb.preU{k};
    uns = find(l < 0 & u > 0);
    x0 = (lb + ub) / 2;
    [~, Z0] = i_forward(ops, x0, reluIdx);
    z0 = Z0{k};
    fx = zeros(numel(l), 1, prec);
    for i = 1:min(4, numel(uns))
        j = uns(i); s = sign(z0(j)); if s == 0, s = 1; end
        fx(j) = s;                              % fix to x0's actual sign -> x0 satisfies it
    end
    fixings{k} = fx;
end

function [Y, Z] = i_forward(ops, X, reluIdx) %#ok<INUSD>
% Exact forward eval of the op list on columns of X (n x m). Y = output (nOut x m); Z{k} = the
% pre-activation at each relu op (for the split-constraint check).
    nOps = numel(ops); co = cell(nOps+1,1); co{1} = X; Z = cell(nOps,1);
    for k = 1:nOps
        op = ops{k};
        if strcmp(op.type,'add')
            a = op.inputs(1)+1; b = op.inputs(2)+1; co{k+1} = co{a} + co{b}; continue;
        end
        x = co{op.src+1};
        switch op.type
            case 'affine',     co{k+1} = op.W*x + op.b(:);
            case 'conv',       co{k+1} = i_conv_fwd(op, x);
            case 'normaffine', co{k+1} = i_bcast(op.scale,op.shape).*x + i_bcast(op.shift,op.shape);
            case 'avgpool',    co{k+1} = i_avgpool_fwd(op, x);
            case 'relu',       Z{k} = x; co{k+1} = max(x,0);
        end
    end
    Y = co{nOps+1};
end

function y = i_conv_fwd(op, x)
    ish=op.inShape; osh=op.outShape; m=size(x,2); bb=reshape(op.b(:),[1 1 osh(3)]);
    X4=dlarray(reshape(x,[ish(1) ish(2) ish(3) m]),'SSCB');
    pad2=[op.pad(1) op.pad(3); op.pad(2) op.pad(4)];
    Y4=dlconv(X4,op.W,bb,'Stride',op.stride,'Padding',pad2,'DilationFactor',op.dil);
    y=reshape(extractdata(Y4),[prod(osh) m]);
end

function y = i_avgpool_fwd(op, x)
    ish=op.inShape; osh=op.outShape; m=size(x,2); kh=op.pool(1); kw=op.pool(2);
    X4=reshape(x,[ish(1) ish(2) ish(3) m]); Y=zeros([osh(1) osh(2) osh(3) m]);
    for oh=1:osh(1)
        for ow=1:osh(2)
            rh=(oh-1)*op.stride(1)+(1:kh); rw=(ow-1)*op.stride(2)+(1:kw);
            Y(oh,ow,:,:)=mean(mean(X4(rh,rw,:,:),1),2);
        end
    end
    y=reshape(Y,[prod(osh) m]);
end

function ok = i_satisfies(Z, fixings, reluIdx, tol)
% Logical (1 x m): each sample satisfies every split constraint (active fix z>=0, inactive z<=0).
    m = size(Z{reluIdx(1)}, 2); ok = true(1,m);
    for r = 1:numel(reluIdx)
        k = reluIdx(r);
        if isempty(fixings) || numel(fixings)<k || isempty(fixings{k}), continue; end
        fx = fixings{k}(:,1); z = Z{k};
        if any(fx==1),  ok = ok & all(z(fx==1,:)  >= -tol, 1); end
        if any(fx==-1), ok = ok & all(z(fx==-1,:) <=  tol, 1); end
    end
end
