
% Test a step Reach with a single input
I = ExamplePoly.randVrep;   % input set
B = outerApprox(I); % interval hull of the input set
lb = B.Internal.lb; % min-vec of x vector
ub = B.Internal.ub; % max-vec of x vector

O1 = ReLU.stepReach(I, 1, lb(1), ub(1), 'exact'); % step reach set corresponding to index = 1

Out = ReLU.reach(I, 'exact');    % reach set of ReLU(I)

fig1 = figure;
I.plot; % plot input set
fig2 = figure;
O1.plot;    % plot output set of ReLU_1 stepReach operation
fig3 = figure; 
Out.plot; % plot output set of ReLU(I)

Out2 = ReLU.reach(I, 'approx'); % over-approximate reach set
fig4 = figure;
Out2.plot;

