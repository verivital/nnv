
I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star


B = I.getBox;
l = B.lb;
u = B.ub;
y_l = logsig(l);
y_u = logsig(u);
dy_l = logsig('dn', l);
dy_u = logsig('dn', u);

S = LogSig.stepLogSig_Split(I, 1, l(1), u(1), y_l(1), y_u(1), dy_l(1), dy_u(1));

figure;
Star.plots(S);

S1 = LogSig.stepLogSig_NoSplit(I, 1, l(1), u(1), y_l(1), y_u(1), dy_l(1), dy_u(1));

figure;
Star.plots(S1);

figure;
Star.plots(S1);
hold on;
Star.plots(S, 'red');


% figure;
% I.plot; % input set
% hold on;
% S.plot; % reach set
% hold on;
% plot(Y(1, :), Y(2, :), '*'); % sampled outputs
% hold on;
% plot(X(1, :), X(2, :), 'ob'); % sampled inputs
% 
