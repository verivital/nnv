lb = [49; 25; 9; 20];
ub = [51; 25.2; 11; 20.2];

B1 = Box(lb, ub);

init_set = B1.toStar();

tic;
B2 = init_set.getBox();
toc;

