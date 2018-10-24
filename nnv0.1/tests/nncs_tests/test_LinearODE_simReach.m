A = [0 1;-5 -2];
x0 = [1; 2];
h = 0.1;
N = 5;

lb = [-1 ;-1];
ub = [1; 1]; 
B = Box(lb, ub);
X0 = B.toStar(); % initial set

R1 = LinearODE.simReachDirect(A, X0, h, N);
R2 = LinearODE.simReachKrylov(A, X0, h, N, 2);
R3 = LinearODE.simReachOde45(A, X0, h, N);

fig = figure;
Star.plots(R1);
fig = figure;
Star.plots(R2);
fig = figure;
Star.plots(R3);