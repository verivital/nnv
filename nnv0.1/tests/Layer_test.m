I = ExamplePoly.randVrep;   % input set
W = [1.5 1; 0 0.5];
b = [0.5; 0.5];
L = Layer(W, b, 'ReLU');
[R, rn, t] = L.reach(I, 'exact');    % exact reach set
[R1, rn1, t1] = L.reach(I, 'approx-oriented-box');  % over-approximate reach set
[R2, rn2, t2] = L.reach(I, 'approx-box'); %over-approximate reach set using box
[R3, rn3, t3] = L.reach(I, 'approx-polyhedron'); %over-approximate reach set using polyhedron

fig = figure; 
I.plot;
title('Input Set');

fig = figure;
R.plot;
title('Exact reach set');

fig = figure;
R1.plot;
title('Over-approximate reach set using oriented box');

fig = figure;
R2.plot;
title('Over-approximate reach set using box');

fig = figure;
R3.plot;
title('Over-approximate reach set using polyhedron');
