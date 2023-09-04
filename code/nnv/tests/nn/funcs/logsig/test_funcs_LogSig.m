% Create random set

I = ExamplePoly.randVrep;   
I.outerApprox;
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b, I.Internal.lb, I.Internal.ub); % input star

B = I.getBox; 
I_zono = B.toZono;

W = [0.5 1; -1 1];
I_zono = I_zono.affineMap(W, []);
I_star = I_zono.toStar;

sample_size=25;

l = B.lb;%these and below are actually not used.
u = B.ub;
y_l = logsig(l);
y_u = logsig(u);
dy_l = logsig('dn', l);
dy_u = logsig('dn', u);


%% test 1: LogSig reach approx star
S = LogSig.reach(I_star, 'approx-star');

%% test 2: LogSig reach approx star split
S = LogSig.reach(I_star, 'approx-star-split');

%% test 3: LogSig reach abs dom
S = LogSig.reach(I_star, 'abs-dom');

%% test 4: LogSig reach approx zono
Z = LogSig.reach(I_zono, 'approx-zono');

%% test 5: LogSig reach zono approx split
S = LogSig.reach_star_approx_split(I_star);

%% test 6: LogSig reach zono approx no spit
S = LogSig.reach_star_approx_no_split(I_star);

%% test 7: LogSig reach zono approx 
Z = LogSig.reach_zono_approx(I_zono);
