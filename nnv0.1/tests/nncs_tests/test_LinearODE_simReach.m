A = [0 1;-5 -2];
x0 = [1; 2];
h = 0.1;
N = 5;

lb = [-1 ;-1];
ub = [1; 1]; 
B = Box(lb, ub);
X0 = B.toStar(); % initial set

R = LinearODE.simReachDirect(A, X0, h, N);
Star.plots(R);