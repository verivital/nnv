I = ExamplePoly.randVrep;   % input set1
I.minVRep();
P = I.V';
V = [0 1; 1 0; 1/sqrt(2) 1/sqrt(2); -1/sqrt(2) 1/sqrt(2)]';    % unit vector

Z1 = Reduction.zonotopeHull(P, V); % specified zonotope generator
Z2 = Reduction.zonotopeHull_Random(P, 20);    % random zonotop generator

fig = figure;
Z1.plot;
hold on;
plot(P(1, :), P(2, :), '*');
hold on;
I.plot;

fig = figure;
Z2.plot;
hold on;
plot(P(1, :), P(2, :), '*');


N = 20;
I_vol = [];
Z1_vol = [];
Z2_vol = [];
I_num_constr = []; % number of constraints of input set
Z1_num_constr = []; % number of constraints of zonotope Z1
Z2_num_constr = []; % number of constraints of zonotope Z2
Z1_vol_inc = [];   % percent of volume increase by over-approximation using box Z1
Z2_vol_inc = [];   % percent of volum increase by over-approximation using zonotope Z2
Z1_num_constr_dec = []; % percent of number of constraints reduced using box Z1
Z2_num_constr_dec = [];

for i=1:N
    
    I = ExamplePoly.randVrep;   % input set
    I.minVRep();
    P = I.V';
    V = [0 1; 1 0]';    % unit vector
    Z1 = Reduction.zonotopeHull(P, V); % specified zonotope generator
    Z2 = Reduction.zonotopeHull_Random(P, 20);    % random zonotop generator
    I_vol = [I_vol I.volume];
    Z1_vol = [Z1_vol Z1.volume];
    Z2_vol = [Z2_vol Z2.volume];
    
    Z1_vol_inc = [Z1_vol_inc (Z1.volume - I.volume)/I.volume];
    Z2_vol_inc = [Z2_vol_inc (Z2.volume - I.volume)/I.volume];
    
    Z1.minHRep();
    Z2.minHRep();
    n0 = size(I.A, 1);
    n1 = size(Z1.A, 1);
    n2 = size(Z2.A, 1);
    I_num_constr = [I_num_constr n0];
    Z1_num_constr = [Z1_num_constr n1];
    Z2_num_constr = [Z2_num_constr n2];
    
    Z1_num_constr_dec = [Z1_num_constr_dec (n0 - n1)/ n0];
    Z2_num_constr_dec = [Z2_num_constr_dec (n0 - n2)/ n0];
end

fig = figure;
x = 1:N;
y1 = Z1_vol_inc;
y2 = Z2_vol_inc;

plot(x, y1, 'g--*');
hold on;
plot(x, y2, 'b--o');
legend('Z1\_vol\_inc', 'Z2\_vol\_inc');
xlabel('N');
ylabel('Percentage of Volume Increased');
title('Percentage of Volume Encreased by Over-Approximation')
saveas(fig, 'ZonotopeHullVolumePerformance.pdf');


fig = figure;
x = 1:N;
y1 = Z1_num_constr_dec;
y2 = Z2_num_constr_dec;

plot(x, y1, 'g--*');
hold on;
plot(x, y2, 'b--o');
legend('Z1\_num\_constr\_dec', 'Z2\_num\_constr\_dec');
xlabel('N');
ylabel('Percentage of number of constraints decreased');
title('Percentage of Constraints decreased by Over-Approximation')
saveas(fig, 'ZonotopeHullConstraintsPerformance.pdf');

fig = figure;
x = 1:N;
y1 = I_num_constr;
y2 = Z1_num_constr;
y3 = Z2_num_constr;

plot(x, y1, 'g--*');
hold on;
plot(x, y2, 'b--o');
hold on;
plot(x, y3, 'r--x');

legend('I', 'Z1', 'Z2');
xlabel('N'); 
ylabel('Number of constraints');
title('Number of constraints');
saveas(fig, 'ZonotopeHullNumberOfConstraints.pdf');
