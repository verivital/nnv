%% test_SoftmaxAttn_softmax - TDD for the exact correlated row-softmax bound
% Mirrors n2v tests/soundness/test_soundness_softmax_attention.py.
% softmax is taken over the KEY axis (dim 2 of an [R x n] logit box).
% Monotonicity: dA_j/dS_j >= 0, dA_j/dS_k <= 0 (k~=j), so the per-coordinate
% extremes are attained at the corners
%   a_ub_j = softmax(s_lo everywhere, s_hi at j)_j
%   a_lb_j = softmax(s_hi everywhere, s_lo at j)_j
% giving the EXACT (tightest) axis-aligned enclosure of softmax over the box.

here = fileparts(mfilename('fullpath'));
nnvroot = fullfile(here, '..', '..', '..');
addpath(genpath(fullfile(nnvroot, 'engine')));
addpath(genpath(fullfile(nnvroot, 'tests', 'soundness')));
rng(20260625);

%% edge cases: equal logits -> uniform 1/n; single key -> exactly 1
s = [2 2 2; -1 -1 -1];               % R=2, n=3, degenerate box (lo==hi)
[al, au] = SoftmaxAttn.correlatedRowSoftmaxBounds(s, s);
assert(max(abs(al(:) - 1/3)) < 1e-9, 'equal logits lower');
assert(max(abs(au(:) - 1/3)) < 1e-9, 'equal logits upper');

[al, au] = SoftmaxAttn.correlatedRowSoftmaxBounds([0.5; -2], [3; 1]);  % R=2, n=1
assert(all(abs(al - 1) < 1e-12) && all(abs(au - 1) < 1e-12), 'single key = 1');

%% soundness + tightness: sampled row-softmax lies in [a_lb,a_ub], bounds attained
regimes = {'pos','neg','wide','mixed'};
shapes  = {[5 5],[17 17],[3 1],[1 6],[4 3]};
for ri = 1:numel(regimes)
  for si = 1:numel(shapes)
    R = shapes{si}(1); n = shapes{si}(2);
    [s_lo, s_hi] = rand_logit_box(R, n, regimes{ri});
    [al, au] = SoftmaxAttn.correlatedRowSoftmaxBounds(s_lo, s_hi);
    assert(all(al(:) >= -1e-12) && all(au(:) <= 1+1e-12), 'within [0,1]');
    assert(all(al(:) <= au(:) + 1e-12), 'lb<=ub');
    for t = 1:400
      S = s_lo + (s_hi - s_lo).*rand(R, n);
      A = row_softmax(S);
      assert(all(A(:) >= al(:) - 1e-9) && all(A(:) <= au(:) + 1e-9), ...
        'softmax UNSOUND regime %s shape %s', regimes{ri}, mat2str(shapes{si}));
    end
    % each row must sum to 1; bound must be consistent with that
    assert(all(sum(al,2) <= 1 + 1e-9) && all(sum(au,2) >= 1 - 1e-9), 'row-sum sandwich');
  end
end
fprintf('PASS correlatedRowSoftmaxBounds soundness+tightness\n');

%% softmaxAttnStar: Star logits -> box Star weights; sampled weights contained
R = 5; n = 5;
Slog = rand_star(R*n, 0.6);
W = SoftmaxAttn.softmaxAttnStar(Slog, [R n], 'estimate');
assert(W.dim == R*n);
for t = 1:200
  a = sample_alpha(Slog); v = Slog.V(:,1) + Slog.V(:,2:end)*a;
  Smat = reshape(v, [R n]);
  A = reshape(row_softmax(Smat), [], 1);
  assert(soundness_test_utils.verify_star_containment(W, A, 1e-6), 'softmaxAttnStar UNSOUND %d', t);
end
fprintf('PASS softmaxAttnStar soundness\n');

% ================= local helpers =================
function A = row_softmax(S)
  m = max(S, [], 2);
  E = exp(S - m);
  A = E ./ sum(E, 2);
end

function [lo, hi] = rand_logit_box(R, n, regime)
  w = 0.3 + 0.8*rand(R, n);
  switch regime
    case 'pos',   lo = 0.2 + rand(R,n);   hi = lo + w;
    case 'neg',   hi = -0.2 - rand(R,n);  lo = hi - w;
    case 'wide',  lo = -3*rand(R,n)-1;    hi = 3*rand(R,n)+1;
    case 'mixed', lo = -rand(R,n)-0.1;    hi = rand(R,n)+0.1;
  end
end

function S = rand_star(d, rad)
  c = randn(d,1);
  S = Star(c - rad, c + rad);
end

function a = sample_alpha(S)
  a = S.predicate_lb + (S.predicate_ub - S.predicate_lb).*rand(S.nVar,1);
  tries = 0;
  while ~isempty(S.C) && any(S.C*a > S.d) && tries < 100
    a = S.predicate_lb + (S.predicate_ub - S.predicate_lb).*rand(S.nVar,1); tries = tries+1;
  end
end
