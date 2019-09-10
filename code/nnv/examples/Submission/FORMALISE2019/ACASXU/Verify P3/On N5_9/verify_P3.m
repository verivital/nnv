load result.mat;
load ACASXU_run2a_5_9_batch_2000.mat;

R = F.outputSet;
n = length(R); 

lb = [1500; -0.06; 3.1; 980; 960];
ub = [1800; 0.06; 3.14; 1200; 1200];

% normalize input
for i=1:5
    lb(i) = (lb(i) - means_for_scaling(i))/range_for_scaling(i);
    ub(i) = (ub(i) - means_for_scaling(i))/range_for_scaling(i);   
end


n=5000;

x1 = (ub(1) - lb(1)).*rand(n, 1) + lb(1);
x2 = (ub(2) - lb(2)).*rand(n, 1) + lb(2);
x3 = (ub(3) - lb(3)).*rand(n, 1) + lb(3);
x4 = (ub(4) - lb(4)).*rand(n, 1) + lb(4);
x5 = (ub(5) - lb(5)).*rand(n, 1) + lb(5);

I = [x1'; x2'; x3'; x4'; x5'];

Y = F.sample(I);

output = Y{1, 7};

% normalize output
n = length(output);
normalized_output = [];
for i=1:n
    out_i = output(:, i);
    x1 = out_i(1) * range_for_scaling(6) + means_for_scaling(6);
    x2 = out_i(2) * range_for_scaling(6) + means_for_scaling(6);
    x3 = out_i(3) * range_for_scaling(6) + means_for_scaling(6);
    x4 = out_i(4) * range_for_scaling(6) + means_for_scaling(6);
    x5 = out_i(5) * range_for_scaling(6) + means_for_scaling(6);
    x = [x1; x2; x3; x4; x5];
    normalized_output = [normalized_output x];
    
end

maps = [1 0 0 0 0; 0 0 0 0 1; 0 0 1 0 0];
output_mapped = maps*normalized_output;

normalized_mat = range_for_scaling(6) * eye(5);
normalized_vec = means_for_scaling(6) * ones(5,1);

R_normalized = [];
for i=1:length(R)
    R_normalized = [R_normalized Reduction.affineMap(R(i), normalized_mat) + normalized_vec];   
end

R1 = [];
for i=1:length(R_normalized)
    R1 = [R1 Reduction.affineMap(R_normalized(i),maps)];   
end

% plot reachable set
fig = figure;
R1.plot;
hold on;
plot3(output_mapped(1, :), output_mapped(2, :), output_mapped(3, :), '*');
xlabel('COC', 'Fontsize', 30);
ylabel('Strong-Right', 'Fontsize', 30);
zlabel('Weak-Right', 'Fontsize', 30);
set(gca, 'Fontsize', 25);

saveas(gcf,'COC_WeakRight_StrongRight.pdf');


% verify safety: COC is not the minimal score 
% unsafe: x1 <= x2, x1 <= x3, x1 <= x4, x1 <= x5

A = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
b = [0; 0; 0; 0];

U = Polyhedron('A', A, 'b', b); % unsafe set

% this property is hold, safe = true == to Reluplex's result
[safe, check_time] = Verifier.checkSafety(R_normalized, U); % verify safety

if safe
    fprintf('\nThe property 3 holds for N5_9, -> safe');
else
    fprintf('\nThe property 3 does not hold for N5_9 -> unsafe');
end 

