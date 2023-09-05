A = [0 1;-5 -2];
x0 = [1; 2];
h = 0.01;
N = 5;

X1 = LinearODE.simDirect(A, x0, h, N);
X2 = LinearODE.simKrylov(A, x0, h, N, 2);
X3 = LinearODE.simOde45(A, x0, h, N);
