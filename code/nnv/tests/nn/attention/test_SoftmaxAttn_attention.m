%% test_SoftmaxAttn_attention - TDD for the single-head and multi-head attention reach
% singleHeadAttn(Q,K,V,scale,[N D],mode): sound reach of softmax(scale*Q*K')*V with
%   Q,K,V matrix-Stars (dim N*D). Q,K used only for bounds (scores); V kept symbolic.
% selfAttentionReach(X,P,mode): full multi-head self-attention with Q/K/V/out
%   projections, head split, concat -- the ViT encoder attention sublayer (no residual).
% Soundness oracle: linprog containment; inputs sampled, true forward must be contained.

here = fileparts(mfilename('fullpath'));
nnvroot = fullfile(here, '..', '..', '..');
addpath(genpath(fullfile(nnvroot, 'engine')));
addpath(genpath(fullfile(nnvroot, 'tests', 'soundness')));
rng(20260626);

%% singleHeadAttn soundness: independent Q,K,V samples -> output contained
for trial = 1:6
  N = randi([2 6]); D = randi([2 16]); scale = 1/sqrt(D);
  Q = rand_box_star(N*D, 0.4); K = rand_box_star(N*D, 0.4); V = rand_box_star(N*D, 0.4);
  O = SoftmaxAttn.singleHeadAttn(Q, K, V, scale, [N D], 'estimate');
  assert(O.dim == N*D);
  for t = 1:120
    q = realize(Q); k = realize(K); v = realize(V);
    Qm = reshape(q,[N D]); Km = reshape(k,[N D]); Vm = reshape(v,[N D]);
    scores = scale * (Qm * Km');
    A = row_softmax(scores);
    Om = reshape(A * Vm, [], 1);
    assert(soundness_test_utils.verify_star_containment(O, Om, 1e-6), ...
      'singleHeadAttn UNSOUND trial %d sample %d (N=%d D=%d)', trial, t, N, D);
  end
end
fprintf('PASS singleHeadAttn soundness\n');

%% selfAttentionReach soundness: small config -> rigorous LP containment
% small N keeps the containment LP cheap; this is the rigorous soundness gate.
cfg = struct('N',4,'E',12,'H',2);
N=cfg.N; E=cfg.E; H=cfg.H; D=E/H; scale=1/sqrt(D);
P = rand_attn_params(E, H, scale);
X = rand_box_star(N*E, 0.05);
AO = SoftmaxAttn.selfAttentionReach(X, N, P, 'estimate');
assert(AO.dim == N*E && AO.nVar >= X.nVar, 'attn out shape/prefix');
for t = 1:60
  x = realize(X);  out = self_attn_eval(reshape(x,[N E]), N, P);
  assert(soundness_test_utils.verify_star_containment(AO, reshape(out,[],1), 1e-5), ...
    'selfAttentionReach UNSOUND (LP) sample %d', t);
end
Y = SoftmaxAttn.prefixAdd(X, AO);
for t = 1:40
  x = realize(X); out = self_attn_eval(reshape(x,[N E]), N, P);
  assert(soundness_test_utils.verify_star_containment(Y, x + reshape(out,[],1), 1e-5), 'residual UNSOUND %d', t);
end
fprintf('PASS selfAttentionReach soundness + residual (LP, small)\n');

%% selfAttentionReach soundness: ViT shapes -> fast estimate-box containment
% For large stars the LP containment is too slow; instead check every true forward
% lies inside the cheap estimate over-approx box. The estimate box OUTER-bounds the
% reachable set, so a sample escaping it is a definite soundness violation.
configs = { struct('N',5,'E',48,'H',3), struct('N',17,'E',48,'H',3) };
for ci = 1:numel(configs)
  cfg = configs{ci};
  N = cfg.N; E = cfg.E; H = cfg.H; D = E/H; scale = 1/sqrt(D);
  P = rand_attn_params(E, H, scale);
  X = rand_box_star(N*E, 0.04);
  AO = SoftmaxAttn.selfAttentionReach(X, N, P, 'estimate');
  assert(AO.dim == N*E, 'attn out dim');
  [albx, aubx] = AO.estimateRanges();
  for t = 1:80
    x = realize(X);  out = reshape(self_attn_eval(reshape(x,[N E]), N, P), [], 1);
    assert(all(out >= albx - 1e-6) && all(out <= aubx + 1e-6), ...
      'selfAttentionReach escapes estimate box (UNSOUND) cfg %d sample %d', ci, t);
  end
end
fprintf('PASS selfAttentionReach soundness (estimate-box, ViT shapes)\n');

% ================= local helpers =================
function A = row_softmax(S)
  m = max(S,[],2); E = exp(S-m); A = E./sum(E,2);
end

function S = rand_box_star(d, rad)
  c = randn(d,1); S = Star(c-rad, c+rad);
end

function v = realize(S)
  a = S.predicate_lb + (S.predicate_ub - S.predicate_lb).*rand(S.nVar,1);
  tries=0;
  while ~isempty(S.C) && any(S.C*a > S.d) && tries<200
    a = S.predicate_lb + (S.predicate_ub - S.predicate_lb).*rand(S.nVar,1); tries=tries+1;
  end
  v = S.V(:,1) + S.V(:,2:end)*a;
end

function P = rand_attn_params(E, H, scale)
  P.E=E; P.H=H; P.D=E/H; P.scale=scale;
  P.Mq=randn(E,E)*0.1; P.bq=randn(E,1)*0.1;
  P.Mk=randn(E,E)*0.1; P.bk=randn(E,1)*0.1;
  P.Mv=randn(E,E)*0.1; P.bv=randn(E,1)*0.1;
  P.Mo=randn(E,E)*0.1; P.bo=randn(E,1)*0.1;
end

function out = self_attn_eval(Xt, N, P)
  % concrete forward of multi-head self-attention; Xt [N,E]
  E=P.E; H=P.H; D=P.D;
  Q = Xt*P.Mq + P.bq';  K = Xt*P.Mk + P.bk';  V = Xt*P.Mv + P.bv';
  O = zeros(N,E);
  for h=1:H
    cols = (h-1)*D + (1:D);
    Qh=Q(:,cols); Kh=K(:,cols); Vh=V(:,cols);
    A = row_softmax(P.scale*(Qh*Kh'));
    O(:,cols) = A*Vh;
  end
  out = O*P.Mo + P.bo';
end
