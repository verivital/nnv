load F.mat;
load ACASXU_run2a_1_1_batch_2000.mat;

R = F.outputSet;
n = length(R); 

R_bounded = [];
R_unbounded = [];

for i=1:n
    fprintf('\nChecking R(%d)', i);
    if ~R(i).isBounded
        fprintf('\nR(%d) is not bounded', i);
        R_unbounded = [R_unbounded R(i)];
    else
        R_bounded = [R_bounded R(i)];
    end
end

B = [];
m = length(R_bounded);
for i=1:m
    B = [B R_bounded(i).outerApprox];
end

lb = [55947.69; -3.14; -3.14; 1145; 0];
ub = [60760; 3.14; 3.14; 1200; 60];

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
    x1 = out_i(1) * range_for_scaling(1) + means_for_scaling(1);
    x2 = out_i(2) * range_for_scaling(2) + means_for_scaling(2);
    x3 = out_i(3) * range_for_scaling(3) + means_for_scaling(3);
    x4 = out_i(4) * range_for_scaling(4) + means_for_scaling(4);
    x5 = out_i(5) * range_for_scaling(5) + means_for_scaling(5);
    x = [x1; x2; x3; x4; x5];
    normalized_output = [normalized_output x];
    
end


maps = [0 0 0 1 0; 0 0 0 0 1];
output_mapped = maps*normalized_output;

R1 = [];
for i=1:length(R_bounded)
    R1 = [R1 R_bounded(i).affineMap(maps)];   
end

fig = figure;
R1.plot;
hold on;
plot(output_mapped(1, :), output_mapped(2, :), 'o');


% find a counter example

n = length(output);

counter_O = [];
counter_I = [];

for i=1:n
    
    [temp, idx] = min(normalized_output(:, i));
    if idx ~=1
        fprintf('\nInput %d produce counter example', i);
        fprintf('\nThis is the output of the neural network:')
        display(normalized_output(:, i));
        fprintf('\nThis is the corresponding input for above output');
        display(I(:, i));
        counter_O = [counter_O normalized_output(:, i)];
        counter_I = [counter_I I(:, i)];
    end
end




