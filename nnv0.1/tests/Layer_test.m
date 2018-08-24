I = ExamplePoly.randVrep;   % input set
W = [1.5 1; 0 0.5];
b = [0.5; 0.5];
L = Layer(W, b, 'ReLU');
[R, rn, t] = L.reach(I, 'exact');    % exact reach set
[R1, rn1, t1] = L.reach(I, 'approx');  % over-approximate reach set
fig0 = figure; 
I.plot;
I1 = I.affineMap(W);
I2 = I1 + b; % Wx + b set, x in I
fig1 = figure;
I2.plot;
fig2 = figure;
R.plot;
fig3 = figure;
R1.plot;
% One can verify that computing an over-approximate reach set is
% usually faster than computing an exact reach set, i.e., t1 < t