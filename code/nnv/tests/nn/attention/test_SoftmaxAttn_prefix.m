%% test_SoftmaxAttn_prefix - TDD for prefix-aligned residual add and head concat
% prefixAdd(x, branch): x's predicates are a PREFIX of branch's (branch was built
%   from x by appending relaxation predicates). Returns x+branch over branch's
%   constraint system -- sound AND tight (keeps residual correlation a block-diagonal
%   Minkowski join would drop), and provenance-safe (never the structural-equality
%   fast path that under-approximates independent-but-identical operands).
% prefixConcat({S_h}, n_old): stack states of stars sharing the first n_old preds
%   into one star with block-diagonal fresh tails -> used to assemble attention heads.

here = fileparts(mfilename('fullpath'));
nnvroot = fullfile(here, '..', '..', '..');
addpath(genpath(fullfile(nnvroot, 'engine')));
addpath(genpath(fullfile(nnvroot, 'tests', 'soundness')));
rng(20260625);

%% prefixAdd: y = x + branch(x), sampled sums contained (residual stream)
% Build x over input preds; branch = avEnvelope(...) which appends fresh preds.
K=4; D=5; M=K;   % residual needs x.dim == branch.dim -> use M=K, D s.t. M*D == K*D? use square
% make x.dim == branch.dim: branch dim = M*D, so set x to dim M*D
Vstar = rand_box_star(K*D, 'mixed');
a_lb = rand(M,K)*0.5; a_ub = a_lb + rand(M,K)*0.5;
branch = SoftmaxAttn.avEnvelopeStar(a_lb, a_ub, Vstar, [K D], 'estimate');  % dim M*D, preds [alpha_old, fresh]
% x: a star over the SAME alpha_old prefix only (dim M*D), built as an affine map of Vstar
Wx = randn(M*D, K*D);
x = Vstar.affineMap(Wx, zeros(M*D,1));     % x.nVar == Vstar.nVar == prefix of branch
assert(x.nVar <= branch.nVar, 'x prefix');
y = SoftmaxAttn.prefixAdd(x, branch);
assert(y.dim == branch.dim && y.nVar == branch.nVar, 'y shape');
for t=1:250
  abeta = sample_alpha(branch);                 % full predicate vector of branch
  aold  = abeta(1:Vstar.nVar);
  vvec = Vstar.V(:,1) + Vstar.V(:,2:end)*aold;  % realized V
  Vm = reshape(vvec,[K D]);
  A = a_lb + (a_ub-a_lb).*rand(M,K);
  br_true = reshape(A*Vm,[],1);                 % a valid branch output for this aold
  x_true  = x.V(:,1) + x.V(:,2:end)*aold;
  y_true  = x_true + br_true;
  assert(soundness_test_utils.verify_star_containment(y, y_true, 1e-6), 'prefixAdd UNSOUND %d', t);
end
fprintf('PASS prefixAdd soundness\n');

%% prefixAdd does NOT under-approximate independent-but-identical operands
% (the structural-equality trap: x in [-1,1], plus an independent -x' in [-1,1];
%  a structural "shared alpha" gate would collapse to {0}; true range is [-2,2]).
% prefixAdd is only valid when one is a genuine prefix of the other; here we
% verify the residual path keeps the full range when branch genuinely extends x.
xb = Star(-1, 1);                       % dim 1, 1 pred
br = Star([-1;-1+0], [1;1]);            % placeholder, rebuild properly below
% branch over [alpha_x, fresh] with state = -fresh (independent of x):
Vb = [0, 0, -1];                        % center 0, prefix gen 0, tail gen -1
Cb = [1 0; -1 0; 0 1; 0 -1];           % box on [alpha_x, fresh]
db = [1;1;1;1];
br = Star(Vb, Cb, db, [-1;-1], [1;1]);
y = SoftmaxAttn.prefixAdd(xb, br);
lo = y.getMin(1,'linprog'); hi = y.getMax(1,'linprog');
assert(lo <= -2+1e-6 && hi >= 2-1e-6, 'prefixAdd must keep [-2,2], got [%g,%g]', lo, hi);
fprintf('PASS prefixAdd keeps independent-operand range\n');

%% prefixConcat: assemble attention heads (slice shared V, av per head, concat)
H=3; K=5; D=16; M=5; n_old_dim = K*D;
Vfull = rand_box_star(H*K*D, 'mixed');   % full value projection
Os = cell(1,H); A_lbs=cell(1,H); A_ubs=cell(1,H);
for h=1:H
  idx = (h-1)*K*D + (1:K*D);
  Vh = slice_star(Vfull, idx);           % head h value sub-star (same preds)
  A_lbs{h} = rand(M,K)*0.5; A_ubs{h} = A_lbs{h} + rand(M,K)*0.5;
  Os{h} = SoftmaxAttn.avEnvelopeStar(A_lbs{h}, A_ubs{h}, Vh, [K D], 'estimate');
end
O = SoftmaxAttn.prefixConcat(Os, Vfull.nVar);
assert(O.dim == H*M*D, 'concat dim');
for t=1:150
  aold = sample_alpha(Vfull);
  vfull = Vfull.V(:,1) + Vfull.V(:,2:end)*aold;
  o_true = [];
  for h=1:H
    idx = (h-1)*K*D + (1:K*D);
    Vm = reshape(vfull(idx),[K D]);
    A = A_lbs{h} + (A_ubs{h}-A_lbs{h}).*rand(M,K);
    o_true = [o_true; reshape(A*Vm,[],1)];
  end
  assert(soundness_test_utils.verify_star_containment(O, o_true, 1e-6), 'prefixConcat UNSOUND %d', t);
end
fprintf('PASS prefixConcat soundness (multi-head assembly)\n');

% ================= local helpers =================
function S = rand_box_star(d, regime)
  switch regime
    case 'pos',   lo = 0.2+rand(d,1);  hi = lo+0.3+rand(d,1);
    case 'neg',   hi = -0.2-rand(d,1); lo = hi-0.3-rand(d,1);
    case 'mixed', lo = -rand(d,1)-0.2; hi = rand(d,1)+0.2;
  end
  S = Star(lo,hi);
end

function Sh = slice_star(S, idx)
  % select state rows idx, keep all predicates/constraints
  Vh = S.V(idx, :);
  Sh = Star(Vh, S.C, S.d, S.predicate_lb, S.predicate_ub);
end

function a = sample_alpha(S)
  a = S.predicate_lb + (S.predicate_ub - S.predicate_lb).*rand(S.nVar,1);
  tries = 0;
  while ~isempty(S.C) && any(S.C*a > S.d) && tries < 200
    a = S.predicate_lb + (S.predicate_ub - S.predicate_lb).*rand(S.nVar,1); tries = tries+1;
  end
end
