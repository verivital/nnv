A = [0 1;-5 -2];
B = [0;3];
C = [0 1];
D = 0;
Ts = 0.1;
sys = LinearODE(A, B, C, D);
sysd = sys.c2d(Ts); % convert from continuous to discrete

[u,t] = gensig('square',4,10,0.1);

x0 = [1; 2];


sys2 = sysd.d2c(); % convert from discrete to continuous

sysd.simulate(u, t, x0); % simulate the system
%sysd.initial(x0, 8); % compute response of the system with initial condition
%sysd.step(8); % step response