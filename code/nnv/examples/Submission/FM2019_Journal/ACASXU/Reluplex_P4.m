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

VT(1,1) = 1291;
VT(1,2) = 1267;
VT(1,3) = 1150;
VT(1,4) = 107;
VT(1,5) = 352;
VT(1,6) = 219;
VT(1,7) = 1;
VT(1,8) = 3;
VT(1,9) = 3;

VT(2,1) = 330;
VT(2,2) = 415;
VT(2,3) = 243;
VT(2,4) = 86;
VT(2,5) = 151;
VT(2,6) = 118;
VT(2,7) = 34;
VT(2,8) = 549;
VT(2,9) = 52;

VT(3,1) = 478;
VT(3,2) = 107;
VT(3,3) = 116;
VT(3,4) = 75;
VT(3,5) = 206;
VT(3,6) = 141;
VT(3,7) = 304;
VT(3,8) = 131;
VT(3,9) = 621;

VT(4,1) = 49;
VT(4,2) = 244;
VT(4,3) = 243;
VT(4,4) = 206;
VT(4,5) = 217;
VT(4,6) = 177;
VT(4,7) = 47;
VT(4,8) = 195;
VT(4,9) = 422;

VT(5,1) = 513;
VT(5,2) = 210;
VT(5,3) = 124;
VT(5,4) = 144;
VT(5,5) = 114;
VT(5,6) = 160;
VT(5,7) = 38;
VT(5,8) = 111;
VT(5,9) = 116;

save Reluplex_P4.mat Res VT;

