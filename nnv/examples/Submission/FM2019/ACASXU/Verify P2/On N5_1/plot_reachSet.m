load outputSet.mat;
load ACASXU_run2a_5_1_batch_2000.mat;

normalized_mat = range_for_scaling(6) * eye(5);
normalized_vec = means_for_scaling(6) * ones(5,1);

% normalize output set

fprintf('\nNormalize output set');

n = length(R1);
R1_scaled = [];
parfor i=1:n
    fprintf('\nScaling %d^th output set',i);
    R1_scaled = [R1_scaled  R1(i).affineMap(normalized_mat, normalized_vec)]; % exact normalized reach set
end
R2_scaled = R2.affineMap(normalized_mat, normalized_vec); % over-approximate normalized reach set using star
R3_scaled = R3.affineMap(normalized_mat, normalized_vec); % over-approximate normalized reach set using zonotope
R4_scaled = R4.affineMap(normalized_mat, normalized_vec); % over-approximate normalized reach set using abstract domain


% output = [x1 = COC; x2 = Weak Left; x3 = Weak Right; x4 = Strong Left; x5 = Strong Right]
maps1 = [1 0 0 0 0; 0 1 0 0 0]; % plot a projection on COC, Weak Left
maps2 = [1 0 0 0 0; 0 0 1 0 0]; % plot a projection on COC, Weak Right
maps3 = [1 0 0 0 0; 0 0 0 1 0]; % plot a projection on COC, Strong Left
maps4 = [1 0 0 0 0; 0 0 0 0 1]; % plot a projection on COC, Strong Right

R11 = [];
R12 = [];
R13 = [];
R14 = [];
parfor i=1:n
    R11 = [R11 R1_scaled(i).affineMap(maps1, [])];
    R12 = [R12 R1_scaled(i).affineMap(maps2, [])];
    R13 = [R13 R1_scaled(i).affineMap(maps3, [])];
    R14 = [R14 R1_scaled(i).affineMap(maps4, [])];
end

R21 = R2_scaled.affineMap(maps1, []);
R22 = R2_scaled.affineMap(maps2, []);
R23 = R2_scaled.affineMap(maps3, []);
R24 = R2_scaled.affineMap(maps4, []);

R31 = R3_scaled.affineMap(maps1, []);
R32 = R3_scaled.affineMap(maps2, []);
R33 = R3_scaled.affineMap(maps3, []);
R34 = R3_scaled.affineMap(maps4, []);

R41 = R4_scaled.affineMap(maps1, []);
R42 = R4_scaled.affineMap(maps2, []);
R43 = R4_scaled.affineMap(maps3, []);
R44 = R4_scaled.affineMap(maps4, []);


% plot reachable set

fprintf('\nPlotting reachable set...');

Layers = [];
n = length(b);
for i=1:n - 1
    bi = cell2mat(b(i));
    Wi = cell2mat(W(i));
    Li = LayerS(Wi, bi, 'poslin');
    Layers = [Layers Li];
end
bn = cell2mat(b(n));
Wn = cell2mat(W(n));
Ln = LayerS(Wn, bn, 'purelin');

Layers = [Layers Ln];
F = FFNNS(Layers);

lb = [1500; -0.06; 3.1; 980; 960];
ub = [1800; 0.06; 3.14; 1200; 1200];

% normalize input
for i=1:5
    lb(i) = (lb(i) - means_for_scaling(i))/range_for_scaling(i);
    ub(i) = (ub(i) - means_for_scaling(i))/range_for_scaling(i);   
end

I = Star(lb, ub);

I1 = I.sample(1000);
Y = F.sample(I1);

for i=1:1000
    Y(:,i) = normalized_mat*Y(:,i) + normalized_vec;
end


fig = figure;
Star.plots(R11);
hold on;
plot(Y(1,1:1000), Y(2,1:1000), 'x');

xlabel('COC', 'Fontsize', 16);
ylabel('Weak-Left', 'Fontsize', 16);
set(gca, 'Fontsize', 16);

fig = figure;
R21.plot;
hold on;
plot(Y(1,1:1000), Y(2,1:1000), 'x');

xlabel('COC', 'Fontsize', 16);
ylabel('Weak-Left', 'Fontsize', 16);
set(gca, 'Fontsize', 16);

