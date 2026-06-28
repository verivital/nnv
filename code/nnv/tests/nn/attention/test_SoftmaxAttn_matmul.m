%% test_SoftmaxAttn_matmul - TDD soundness tests for the sound set@set bilinear matmul
% Mirrors n2v tests/soundness/test_soundness_bilinear_matmul.py:
%   - Rump midpoint-radius interval matmul encloses every A@B with A in [al,au], B in [bl,bu]
%   - operands sampled INDEPENDENTLY across sign regimes {pos,neg,mixed} (worst case for
%     concretize, which drops cross-correlation)
%   - bilinearMatMulStar: Star in -> box Star out, every sampled product contained (LP oracle)
%   - negative scale handled sign-safe
%
% Soundness oracle: soundness_test_utils.verify_star_containment (linprog feasibility LP).
% NOTE: tol band here is 0 (exact) for interval bounds, and the util's default band for stars.

% ---- shared setup (runs before each %% section) ----
here = fileparts(mfilename('fullpath'));
nnvroot = fullfile(here, '..', '..', '..');           % code/nnv
addpath(genpath(fullfile(nnvroot, 'engine')));
addpath(genpath(fullfile(nnvroot, 'tests', 'soundness')));
rng(20260625);

%% intervalMatMul encloses all sampled products (sign regimes, multiple shapes)
shapes = {[3 4 5], [5 16 5], [5 5 16], [2 2 4], [1 8 1]};  % [m k n]
regimes = {'pos','neg','mixed'};
for si = 1:numel(shapes)
  s = shapes{si}; m=s(1); k=s(2); n=s(3);
  for ri = 1:numel(regimes)
    [al,au] = rand_box(m,k,regimes{ri});
    [bl,bu] = rand_box(k,n,regimes{ri});
    [cl,cu] = SoftmaxAttn.intervalMatMul(al,au,bl,bu);
    assert(isequal(size(cl),[m n]) && isequal(size(cu),[m n]), 'shape');
    assert(all(cl(:) <= cu(:) + 1e-12), 'lb<=ub');
    % brute-force: sample A,B independently, product must lie in [cl,cu]
    for t = 1:300
      A = al + (au-al).*rand(m,k);
      B = bl + (bu-bl).*rand(k,n);
      C = A*B;
      assert(all(C(:) >= cl(:) - 1e-9) && all(C(:) <= cu(:) + 1e-9), ...
        'interval matmul UNSOUND shape %s regime %s', mat2str(s), regimes{ri});
    end
    % corner check: all-low and all-high operands are inside
    assert(all(vec(al*bl) >= cl(:)-1e-9 & vec(al*bl) <= cu(:)+1e-9));
    assert(all(vec(au*bu) >= cl(:)-1e-9 & vec(au*bu) <= cu(:)+1e-9));
  end
end
fprintf('PASS intervalMatMul soundness\n');

%% bilinearMatMulStar: Star@Star -> box Star, sampled products contained (estimate mode)
m=3; k=4; n=5;
X = rand_star(m*k, 0.5);
Y = rand_star(k*n, 0.5);
S = SoftmaxAttn.bilinearMatMulStar(X, Y, [m k], [k n], 1.0, 'estimate');
assert(S.dim == m*n, 'out dim');
for t = 1:200
  ax = sample_alpha(X); xvec = X.V(:,1) + X.V(:,2:end)*ax;
  ay = sample_alpha(Y); yvec = Y.V(:,1) + Y.V(:,2:end)*ay;
  Xm = reshape(xvec,[m k]); Ym = reshape(yvec,[k n]);
  prod = reshape(Xm*Ym,[],1);
  ok = soundness_test_utils.verify_star_containment(S, prod, 1e-6);
  assert(ok, 'bilinearMatMulStar UNSOUND sample %d', t);
end
fprintf('PASS bilinearMatMulStar soundness (estimate)\n');

%% bilinearMatMulStar with negative scale stays sound
m=2;k=3;n=2; scale=-0.25;
X = rand_star(m*k, 0.7); Y = rand_star(k*n, 0.7);
S = SoftmaxAttn.bilinearMatMulStar(X, Y, [m k], [k n], scale, 'estimate');
for t=1:200
  ax=sample_alpha(X); xvec=X.V(:,1)+X.V(:,2:end)*ax;
  ay=sample_alpha(Y); yvec=Y.V(:,1)+Y.V(:,2:end)*ay;
  prod = reshape(scale*(reshape(xvec,[m k])*reshape(yvec,[k n])),[],1);
  assert(soundness_test_utils.verify_star_containment(S, prod, 1e-6), 'neg-scale UNSOUND %d', t);
end
fprintf('PASS bilinearMatMulStar negative scale\n');

%% lp mode is at least as tight as estimate and still sound
m=3;k=3;n=3;
X = rand_star(m*k, 0.4); Y = rand_star(k*n, 0.4);
Se = SoftmaxAttn.bilinearMatMulStar(X, Y, [m k], [k n], 1.0, 'estimate');
Sl = SoftmaxAttn.bilinearMatMulStar(X, Y, [m k], [k n], 1.0, 'lp');
[le,ue]=Se.getRanges(); [ll,ul]=Sl.getRanges();
assert(all(ll >= le - 1e-7) && all(ul <= ue + 1e-7), 'lp should be tighter-or-equal');
for t=1:100
  ax=sample_alpha(X); xvec=X.V(:,1)+X.V(:,2:end)*ax;
  ay=sample_alpha(Y); yvec=Y.V(:,1)+Y.V(:,2:end)*ay;
  prod = reshape(reshape(xvec,[m k])*reshape(yvec,[k n]),[],1);
  assert(soundness_test_utils.verify_star_containment(Sl, prod, 1e-6), 'lp UNSOUND %d', t);
end
fprintf('PASS bilinearMatMulStar lp tighter+sound\n');

% ================= local helpers =================
function v = vec(M), v = M(:); end

function [lo,hi] = rand_box(r,c,regime)
  w = 0.5 + rand(r,c);          % positive width
  switch regime
    case 'pos',   lo = 0.1 + rand(r,c);        hi = lo + w;
    case 'neg',   hi = -0.1 - rand(r,c);       lo = hi - w;
    case 'mixed', lo = -rand(r,c)-0.1;         hi = rand(r,c)+0.1;
  end
end

function S = rand_star(d, rad)
  % box star [c-r, c+r] with random center, as a Star with explicit predicate box
  c = randn(d,1);
  lo = c - rad; hi = c + rad;
  S = Star(lo, hi);
end

function a = sample_alpha(S)
  a = S.predicate_lb + (S.predicate_ub - S.predicate_lb).*rand(S.nVar,1);
  tries=0;
  while ~isempty(S.C) && any(S.C*a > S.d) && tries<100
    a = S.predicate_lb + (S.predicate_ub - S.predicate_lb).*rand(S.nVar,1); tries=tries+1;
  end
end
