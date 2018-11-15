lb1 = [1; 2; 3; 1];
ub1 = [1.5; 2.5; 3.5; 1.5];

lb2 = [0; 0.5; 1; 1.5];
ub2 = [1; 1; 2; 2];

B1 = Box(lb1, ub1);
B2 = Box(lb2, ub2);

Boxes = [B1 B2];

figure;
PLOT.plotBoxes_2D(Boxes, 1, 2, 'red');

figure;
PLOT.plotBoxes_2D_noFill(Boxes, 1, 2, 'red');

figure;
PLOT.plotBoxes_3D(Boxes, 1, 2, 3, 'red');