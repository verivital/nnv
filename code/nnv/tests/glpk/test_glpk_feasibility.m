
A = [[1, 0]];
b = [0];

c =[0, 0];

glpk(c, A, b);

lb = [];
ub = [];
ctype= [];
vartype = [];
sense = [];
param(1).msglev = 3;
param(1).presol = 0;

[xopt,fopt,status,extra] = glpk(c,A,b,lb, ub, ctype, vartype, sense, param);
status