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

figure;
Star.plots(R1);
figure;
Star.plots(R2);
fig = figure;
Star.plots(R3);


B = [0;3];
C = [0 1];
D = 0;
% set of control input: -0.5 <= u <= 1
lb1 = [-1];
ub1 = [0.5];
Bu = Box(lb1, ub1);
U = Bu.toStar();

sys = LinearODE(A, B, C, D);

R11 = sys.simReach('direct', X0, U, h, N, []); 
R21 = sys.simReach('ode45', X0, U, h, N, []);
R31 = sys.simReach('krylov', X0, U, h, N, 2);
figure;
Star.plots(R11);
figure;
Star.plots(R21);
figure;
Star.plots(R31);