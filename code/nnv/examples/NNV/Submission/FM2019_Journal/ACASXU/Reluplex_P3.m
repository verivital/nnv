VT = zeros(5, 9);
Res = cell(5,9);

Res{1,1} = 'UNSAT';
Res{1,2} = 'UNSAT';
Res{1,3} = 'UNSAT';
Res{1,4} = 'UNSAT';
Res{1,5} = 'UNSAT';
Res{1,6} = 'UNSAT';
Res{1,7} = 'SAT';
Res{1,8} = 'SAT';
Res{1,9} = 'SAT';

Res{2,1} = 'UNSAT';
Res{2,2} = 'UNSAT';
Res{2,3} = 'UNSAT';
Res{2,4} = 'UNSAT';
Res{2,5} = 'UNSAT';
Res{2,6} = 'UNSAT';
Res{2,7} = 'UNSAT';
Res{2,8} = 'UNSAT';
Res{2,9} = 'UNSAT';

Res{3,1} = 'UNSAT';
Res{3,2} = 'UNSAT';
Res{3,3} = 'UNSAT';
Res{3,4} = 'UNSAT';
Res{3,5} = 'UNSAT';
Res{3,6} = 'UNSAT';
Res{3,7} = 'UNSAT';
Res{3,8} = 'UNSAT';
Res{3,9} = 'UNSAT';

Res{4,1} = 'UNSAT';
Res{4,2} = 'UNSAT';
Res{4,3} = 'UNSAT';
Res{4,4} = 'UNSAT';
Res{4,5} = 'UNSAT';
Res{4,6} = 'UNSAT';
Res{4,7} = 'UNSAT';
Res{4,8} = 'UNSAT';
Res{4,9} = 'UNSAT';

Res{5,1} = 'UNSAT';
Res{5,2} = 'UNSAT';
Res{5,3} = 'UNSAT';
Res{5,4} = 'UNSAT';
Res{5,5} = 'UNSAT';
Res{5,6} = 'UNSAT';
Res{5,7} = 'UNSAT';
Res{5,8} = 'UNSAT';
Res{5,9} = 'UNSAT';

VT(1,1) = 6156;
VT(1,2) = 4942;
VT(1,3) = 1134;
VT(1,4) = 528;
VT(1,5) = 317;
VT(1,6) = 64;
VT(1,7) = 1;
VT(1,8) = 3;
VT(1,9) = 2;

VT(2,1) = 1208;
VT(2,2) = 653;
VT(2,3) = 1043;
VT(2,4) = 41;
VT(2,5) = 235;
VT(2,6) = 79;
VT(2,7) = 116;
VT(2,8) = 83;
VT(2,9) = 27;

VT(3,1) = 177;
VT(3,2) = 1561;
VT(3,3) = 1105;
VT(3,4) = 209;
VT(3,5) = 83;
VT(3,6) = 255;
VT(3,7) = 35;
VT(3,8) = 179;
VT(3,9) = 118;

VT(4,1) = 201;
VT(4,2) = 2882;
VT(4,3) = 1767;
VT(4,4) = 86;
VT(4,5) = 35;
VT(4,6) = 300;
VT(4,7) = 126;
VT(4,8) = 142;
VT(4,9) = 143;

VT(5,1) = 1131;
VT(5,2) = 151;
VT(5,3) = 341;
VT(5,4) = 62;
VT(5,5) = 86;
VT(5,6) = 283;
VT(5,7) = 35;
VT(5,8) = 310;
VT(5,9) = 19;

save Reluplex_P3.mat Res VT;

