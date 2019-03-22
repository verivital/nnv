 load mnist_model.mat;
load ones.mat;

Layers = [];
for i=1:6
    W = mnist_model.W{1, i};
    b = mnist_model.b{1,i}';
    L = LayerS(W, b, 'poslin');
    Layers = [Layers L];
end

F = FFNNS(Layers);

% we perform the test on 100 images of each digit
n_inputs = size(ones, 1); % number of digit 1's images in the test: 100


b = 0.0005; % bound of disturbance
db = 0.00001; % increase bound step
k_max = 100; % maximum number of iteration

input_vec = ones(1, :)';
lb = input_vec;
ub = input_vec;
N = length(lb);

t = tic;

k = 0;
bmax = 0;
while(k < k_max)
    
    for i=1:N
        if lb(i) - b > 0
            lb(i) = lb(i) - b;
        else
            lb(i) = 0;
        end
        if ub(i) + b < 1
            ub(i) = ub(i) + b;
        else
            ub(i) = 1;
        end  
    end

    I = Box(lb, ub);

    [R, ~] = F.reach(I.toZono, 'approx-zono');
    B = R.getBox;
    
    if B.lb > 0.5 && B.ub < 1.5
        bmax = b;
        fprintf('\nk = %d -> lb = %.f, ub = %.f, bmax = %.5f', k, B.lb, B.ub, bmax);
        display(B.lb);
        display(B.ub);        
        b = b + db;
        k = k + 1;
    else
        break;
    end
    
end

search_time_zono = toc(t);

% output set value with bmax
input_vec = ones(1, :)';
lb = input_vec;
ub = input_vec;
N = length(lb);
for i=1:N
    if lb(i) - bmax > 0
        lb(i) = lb(i) - bmax;
    else
        lb(i) = 0;
    end
    if ub(i) + bmax < 1
        ub(i) = ub(i) + bmax;
    else
        ub(i) = 1;
    end  
end

I = Box(lb, ub);

[R, ~] = F.reach(I.toZono, 'approx-zono');
B = R.getBox;







