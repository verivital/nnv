%% test_SoftmaxAttn_avenvelope - TDD for the sound symbolic A*V envelope
% Mirrors n2v tests/soundness/test_soundness_av_envelope.py.
% O = A @ V, with A in [a_lb,a_ub] (a_lb>=0, the softmax weights, concretized to a
% box) and V a SYMBOLIC Star. The envelope is affine in V's predicates only, so O
% stays correlated with the input. Sign-aware per-term McCormick (the V<=0 and mixed
% branches are the historically-wrong ones -> hammered here across V regimes).
% A and V are sampled INDEPENDENTLY; every true product must be contained (LP oracle).

here = fileparts(mfilename('fullpath'));
nnvroot = fullfile(here, '..', '..', '..');
addpath(genpath(fullfile(nnvroot, 'engine')));
addpath(genpath(fullfile(nnvroot, 'tests', 'soundness')));
rng(20260625);

%% soundness across V sign regimes (pos/neg/mixed) and shapes; A,V independent
regimes = {'pos','neg','mixed'};
shapes  = {[3 4 5],[5 5 16],[17 17 16],[2 3 2]};   % [M K D]
for ri = 1:numel(regimes)
  for si = 1:numel(shapes)
    M = shapes{si}(1); K = shapes{si}(2); D = shapes{si}(3);
    % attention weights box A in [a_lb,a_ub], a_lb>=0
    a_lb = rand(M,K)*0.5;
    a_ub = a_lb + rand(M,K)*0.5;
    Vstar = rand_value_star(K, D, regimes{ri});
    O = SoftmaxAttn.avEnvelopeStar(a_lb, a_ub, Vstar, [K D], 'estimate');
    assert(O.dim == M*D, 'O dim');
    % O predicates >= V's predicates (V is a prefix)
    assert(O.nVar == Vstar.nVar + M*D, 'fresh predicate count');
    for t = 1:250
      av = sample_alpha(Vstar);
      vvec = Vstar.V(:,1) + Vstar.V(:,2:end)*av;     % realized V
      Vm = reshape(vvec, [K D]);
      A = a_lb + (a_ub - a_lb).*rand(M,K);           % independent A in box
      Om = A * Vm;                                    % true product [M D]
      ovec = reshape(Om, [], 1);
      assert(soundness_test_utils.verify_star_containment(O, ovec, 1e-6), ...
        'avEnvelope UNSOUND regime %s shape %s sample %d', regimes{ri}, mat2str(shapes{si}), t);
    end
  end
end
fprintf('PASS avEnvelopeStar soundness (pos/neg/mixed, independent A,V)\n');

%% degenerate V (point value) and tight A still sound
M=4;K=3;D=2;
c = randn(K*D,1);
Vstar = Star(c, c);              % zero-width V (point)
a_lb = rand(M,K); a_ub = a_lb;   % point A
O = SoftmaxAttn.avEnvelopeStar(a_lb, a_ub, Vstar, [K D], 'estimate');
Vm = reshape(c,[K D]); Om = reshape(a_lb*Vm,[],1);
assert(soundness_test_utils.verify_star_containment(O, Om, 1e-5), 'degenerate UNSOUND');
fprintf('PASS avEnvelopeStar degenerate\n');

%% negative a_lb must raise (the floor>=0 contract)
M=2;K=2;D=2;
Vstar = rand_value_star(K,D,'mixed');
threw = false;
try
  SoftmaxAttn.avEnvelopeStar(-0.1*ones(M,K), ones(M,K), Vstar, [K D], 'estimate');
catch ME
  threw = contains(ME.identifier,'SoftmaxAttn') || contains(ME.message,'a_lb');
end
assert(threw, 'negative a_lb should raise');
fprintf('PASS avEnvelopeStar rejects negative a_lb\n');

% ================= local helpers =================
function S = rand_value_star(K, D, regime)
  % random box Star for V (dim K*D), forced into a sign regime
  d = K*D;
  switch regime
    case 'pos',   lo = 0.2 + rand(d,1);    hi = lo + 0.3 + rand(d,1);
    case 'neg',   hi = -0.2 - rand(d,1);   lo = hi - 0.3 - rand(d,1);
    case 'mixed', lo = -rand(d,1)-0.2;     hi = rand(d,1)+0.2;
  end
  S = Star(lo, hi);
end

function a = sample_alpha(S)
  a = S.predicate_lb + (S.predicate_ub - S.predicate_lb).*rand(S.nVar,1);
  tries = 0;
  while ~isempty(S.C) && any(S.C*a > S.d) && tries < 100
    a = S.predicate_lb + (S.predicate_ub - S.predicate_lb).*rand(S.nVar,1); tries = tries+1;
  end
end
