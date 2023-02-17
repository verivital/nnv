load net256x6.mat;
load inputSet.mat;
N = 25;

% verify the network with eps = 0.02
[r1, rb1, cE1, cands1, vt1] = net.evaluateRBN(S_eps_002(1:N), labels(1:N), 'approx-star', 6);

% % verify the network with eps = 0.05
% [r2, rb2, cE2, cands2, vt2] = net.evaluateRBN(S_eps_005(1:N), labels(1:N), 'approx-star', 6);

% % buid table 

epsilon = [0.02];
verify_time = [sum(vt1)];
safe = [sum(rb1==1)];
unsafe = [sum(rb1 == 0)];
unknown = [sum(rb1 == 2)];


% epsilon = [0.02; 0.05];
% verify_time = [sum(vt1); sum(vt2)];
% safe = [sum(rb1==1); sum(rb2 == 1)];
% unsafe = [sum(rb1 == 0); sum(rb2 == 0)];
% unknown = [sum(rb1 == 2); sum(rb2 == 2)];

T = table(epsilon, safe, unsafe, unknown, verify_time)

save verify_net256x6.mat T r1 rb1 cE1 cands1 vt1 %r2 rb2 cE2 cands2 vt2; 