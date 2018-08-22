I = ExamplePoly.randVrep;   % input set
W = [1.5 1; 0 0.5; 1 -1];
b = [0.5; 0.5; 0];
L = Layer(W, b, 'ReLU');
[R, rn, t] = L.reach(I, 'exact');    % exact reach set
[R1, rn1, t1] = L.reach(I, 'approx');  % over-approximate reach set
fig1 = figure;
R.plot;
fig2 = figure;
R1.plot;
% One can verify that computing an over-approximate reach set is
% usually faster than computing an exact reach set, i.e., t1 < t