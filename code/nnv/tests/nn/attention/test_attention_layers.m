%% test_attention_layers - TDD for sound reach in the NNV transformer layers
% DynamicMatmulLayer.reach: sound set@set matmul (was an unconditional error).
% ScaledDotProductAttentionLayer.reach: sound MULTI-TOKEN attention (was box-lift +
%   single-token V-passthrough that ERRORED on multi-token). All star-family solving
%   methods must return a sound set (never error, never an unsound bound).

here = fileparts(mfilename('fullpath'));
nnvroot = fullfile(here, '..', '..', '..');
addpath(genpath(fullfile(nnvroot, 'engine')));
addpath(genpath(fullfile(nnvroot, 'tests', 'soundness')));
rng(20260626);

%% DynamicMatmulLayer sound reach (set@set), every method, soundness contained
m=3; k=4; n=5;
L = DynamicMatmulLayer('dm');
L.LeftShape = [m k]; L.RightShape = [k n];
A = rand_star(m*k, 0.5); B = rand_star(k*n, 0.5);
methods = {'approx-star','exact-star','abs-dom','relax-star-area','approx-zono'};
for mi = 1:numel(methods)
  S = L.reach({A,B}, methods{mi});
  if isa(S,'Zono'), [slb,sub]=S.getBounds(); contain=@(p) all(p>=slb-1e-6 & p<=sub+1e-6);
  else, contain=@(p) soundness_test_utils.verify_star_containment(S, p, 1e-6); end
  for t=1:60
    a=realize(A); b=realize(B);
    prod = reshape(reshape(a,[m k])*reshape(b,[k n]),[],1);
    assert(contain(prod), 'DynamicMatmul UNSOUND method %s sample %d', methods{mi}, t);
  end
end
fprintf('PASS DynamicMatmulLayer sound reach (all methods)\n');

%% DynamicMatmulLayer without shapes still refuses (sound-by-refusal)
L2 = DynamicMatmulLayer('dm2');
threw=false;
try, L2.reach({A,B},'approx-star'); catch, threw=true; end
assert(threw, 'must refuse without operand shapes');
fprintf('PASS DynamicMatmulLayer refuses unknown shapes\n');

%% ScaledDotProductAttentionLayer multi-token sound reach, all methods
N=5; D=8; scale=1/sqrt(D);
L = ScaledDotProductAttentionLayer('sdpa', D, D, N);   % QueryDim=KeyDim=D, ValueDim=D, SeqLength=N
Q = rand_star(N*D, 0.4); K = rand_star(N*D, 0.4); V = rand_star(N*D, 0.4);
methods = {'approx-star','exact-star','abs-dom','relax-star-area'};
for mi=1:numel(methods)
  S = L.reach(Q, K, V, methods{mi});
  assert(isa(S,'Star') && S.dim==N*D, 'attn out shape (%s)', methods{mi});
  for t=1:50
    q=realize(Q); k_=realize(K); v=realize(V);
    Qm=reshape(q,[N D]); Km=reshape(k_,[N D]); Vm=reshape(v,[N D]);
    A = row_softmax(scale*(Qm*Km'));
    Om = reshape(A*Vm,[],1);
    assert(soundness_test_utils.verify_star_containment(S, Om, 1e-6), ...
      'SDPA multi-token UNSOUND method %s sample %d', methods{mi}, t);
  end
end
fprintf('PASS ScaledDotProductAttentionLayer multi-token sound (all methods)\n');

%% SDPA single-token contract still exact (V-passthrough), no regression
N1=1; D1=6;
L = ScaledDotProductAttentionLayer('sdpa1', D1, D1, N1);
Q=rand_star(D1,0.3); K=rand_star(D1,0.3); V=rand_star(D1,0.3);
S = L.reach(Q,K,V,'approx-star');
for t=1:50
  v=realize(V);                         % single token: attn = V exactly
  assert(soundness_test_utils.verify_star_containment(S, v, 1e-6), 'single-token regressed %d', t);
end
fprintf('PASS ScaledDotProductAttentionLayer single-token exact\n');

% ---- helpers ----
function A = row_softmax(S), m=max(S,[],2); e=exp(S-m); A=e./sum(e,2); end
function S = rand_star(d, rad), c=randn(d,1); S=Star(c-rad,c+rad); end
function v = realize(S)
  a=S.predicate_lb+(S.predicate_ub-S.predicate_lb).*rand(S.nVar,1); tries=0;
  while ~isempty(S.C) && any(S.C*a>S.d) && tries<200
    a=S.predicate_lb+(S.predicate_ub-S.predicate_lb).*rand(S.nVar,1); tries=tries+1;
  end
  v=S.V(:,1)+S.V(:,2:end)*a;
end
