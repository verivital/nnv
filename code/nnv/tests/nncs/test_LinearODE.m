A = [0 1;-5 -2];
B = [0;3];
C = [0 1];
D = 0;
sys = LinearODE(A, B, C, D);
[u,t] = gensig('square',4,10,0.1);

x0 = [1; 2];

sys.simulate(u, t, x0);
%sys.initial(x0, 8);
%sys.step(8);