W1 = [1 -1; 0.5 2; -1 1];
b1 = [-1; 0.5; 0];

W2 = [-2 1 1; 0.5 1 1];
b2 = [-0.5; -0.5];

W3 = [2 -1; 0 1];
b3 = [1; 0];

L1 = LayerS(W1, b1, 'poslin'); % construct first layer
L2 = LayerS(W2, b2, 'poslin');   % construct second layer

lb = -rand(2, 1); % lower-bound vector of input set
ub = lb + [0.5; 0.5];   % upper-bound vector of input set

I = Star(lb, ub); % construct input set
I1 = I.getZono;

X1 = I.affineMap(W1, b1);
BX1 = X1.getBox;
BX11 = X1.getBoxFast;
Y1 = PosLin.reach_star_approx_fast(X1);
B1 = Y1.getBoxFast;
BY1 = Y1.getBox;

XZ1 = I1.affineMap(W1, b1);
Z1 = PosLin.reach_zono_approx(XZ1);
BZ1 = Z1.getBox;

X2 = Y1.affineMap(W2, b2);
BX2 = X2.getBox;
BX22 = X2.getBoxFast;
Y2 = PosLin.reach_star_approx_fast(X2);
B2 = Y2.getBoxFast;

XZ2 = Z1.affineMap(W2, b2);
BXZ2 = XZ2.getBox;
Z2 = PosLin.reach_zono_approx(XZ2);
BZ2 = Z2.getBox;


X3 = Y2.affineMap(W3, b3);
BX3 = X3.getBox;
Y3 = PosLin.reach_star_approx_fast(X3);

XZ3 = Z2.affineMap(W3, b3);
BXZ3 = XZ3.getBox;
Z3 = PosLin.reach_zono_approx(XZ3);



figure;
Z2.plot;
hold on;
Y2.plot;

figure;
Z3.plot;
hold on;
Y3.plot;

