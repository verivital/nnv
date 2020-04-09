load ACASXU_run2a_4_1_batch_2000.mat;
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

% Input Constraints
% 1500 <= i1(\rho) <= 1800,
% -0.06 <= i2 (\theta) <= 0.06,
% \shi = 0
% 1000 <= i4 (\v_own) <= 1200, 
% 700 <= i5 (\v_in) <= 800

lb = [1500; -0.06; 0; 1000; 700];
ub = [1800; 0.06; 0; 1200; 800];

% normalize input
for i=1:5
    lb(i) = (lb(i) - means_for_scaling(i))/range_for_scaling(i);
    ub(i) = (ub(i) - means_for_scaling(i))/range_for_scaling(i);   
end

I = Star(lb, ub);

numCores = 6;

[R1, ~] = F.reach(I, 'exact-star', numCores); % exact reach set using polyhdedron
F.print('F_exact_star.info'); % print all information to a file

[R2, ~] = F.reach(I, 'approx-star'); % approximate reach set using star
F.print('F_approx_star.info'); % print all information to a file

[R3, ~] = F.reach(I.getZono, 'approx-zono'); % approximate reach set using zonotope
F.print('F_approx_zono.info'); % print all information to a file

[R4, ~] = F.reach(I, 'abs-dom'); % approximate reach set using abstract domain
F.print('F_abs_dom.info'); % print all information to a file

save outputSet.mat R1 R2 R3 R4;

normalized_mat = range_for_scaling(6) * eye(5);
normalized_vec = means_for_scaling(6) * ones(5,1);

% normalize output set

fprintf('\nNormalize output set');

n = length(R1);
R1_scaled = [];
for i=1:n
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
for i=1:n
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

fig = figure;
subplot(1, 4, 1);
% R41.plot;
% hold on;
R31.plot;
hold on;
R21.plot;
hold on;
Star.plots(R11);

xlabel('COC', 'Fontsize', 16);
ylabel('Weak-Left', 'Fontsize', 16);
set(gca, 'Fontsize', 16);

subplot(1, 4, 2);
% R42.plot;
% hold on;
R32.plot;
hold on;
R22.plot;
hold on;
Star.plots(R12);
xlabel('COC', 'Fontsize', 16);
ylabel('Weak-Right', 'Fontsize', 16);
set(gca, 'Fontsize', 16);

subplot(1, 4, 3);
% R43.plot;
% hold on;
R33.plot;
hold on;
R23.plot;
hold on;
Star.plots(R13);
xlabel('COC', 'Fontsize', 16);
ylabel('Strong-Left', 'Fontsize', 16);
set(gca, 'Fontsize', 16);

subplot(1, 4, 4);
R34.plot;
hold on;
% R44.plot;
% hold on;
R24.plot;
hold on;
Star.plots(R14);
xlabel('COC', 'Fontsize', 16);
ylabel('Strong-Right', 'Fontsize', 16);
set(gca, 'Fontsize', 16);


set(gca, 'Fontsize', 16);
saveas(gcf, 'reachSet_P4_on_N4_1.pdf');
