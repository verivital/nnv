% This script performs rechability analysis for Emergency Braking System
% author: Dung Tran
% date: 6/3/2019

% normalization 
norm_mat = [1/250 0 0; 0 3.6/120 0; 0 0 1/20]; % normalized matrix

% reinforcement learning controller
load controller.mat; % controller with saturation activation function at the output
rl_layer1 = LayerS(W{1, 1}, b{1, 1}', 'poslin');
rl_layer2 = LayerS(W{1, 2}, b{1, 2}', 'poslin');
rl_layer3 = LayerS(W{1, 3}, b{1, 3}', 'satlin');
rl_controller = FFNNS([rl_layer1 rl_layer2 rl_layer3]); 

% transformation 
load transform.mat;
tf_layer1 = LayerS(W{1, 1}, b{1, 1}', 'poslin');
tf_layer2 = LayerS(W{1, 2}, b{1, 2}', 'poslin');
tf_layer3 = LayerS(W{1, 3}, b{1, 3}', 'purelin');
transformer = FFNNS([tf_layer1 tf_layer2 tf_layer3]);

% control signal scale
scale_mat = [-15*120/3.6 15*10];

% plant matrices
A = [1 -1/15 0; 0 1 0; 0 0 0];
B = [0;1/15;1];
C = [1 0 0;0 1 0; 0 0 1];
D = [0;0;0];
Ts = 0.1;
plant = DLinearODE(A, B, C, D, Ts);


% computation steps
% step 0: initial state of plant
% step 1: normalizing state
% step 2: compute brake output of rl_controller
% step 3: compose brake output and normalized speed
% step 4: compute transformer output
% step 5: get control output
% step 6: scale control output
% step 7: compute reachable set of the plant with the new control input and
% initial state
% step 8: go back to step 1, ...



v0_max = 25.5;
v0_min = 25.2;
d0_max = 55;
a0_min = 0;
a0_max = 0;

d0 = 1; 

l_max = 2; % maximum number of searching 
l = 1;
critical_v0 = [];
while(l < l_max)
    d0_max = d0_max - d0;
    d0_min = d0_max - d0;
    critical_d0 = d0_max;
    inv_tau = 25/(2*v0_max); % worst case full braking time
    num_steps = ceil(15/(inv_tau)) + 20; % minimum number of time steps require for analysis
    lb = [d0_min; v0_min; a0_min]; 
    ub = [d0_max; v0_max; a0_max];
    init_set = Star(lb, ub);
    S = get_reach_set(norm_mat, transformer, rl_controller, scale_mat, A, B, init_set, num_steps);
    inv_TTC_max = get_inv_TTC_max(S);
    if inv_TTC_max >= inv_tau
        fprintf('\nsafety critical intial condition for d_0 = %.5f is found at l = %d', critical_d0, l);
        break;
    end 
end

function S = get_reach_set(norm_mat, transformer, rl_controller, scale_mat, A, B, init_set, num_steps)
    N = num_steps; % number of control steps
    X0 = init_set; % step 0: initial state of plant
    S = X0;
    for i=1:N
        % step 1: normalizing state
        norm_X = X0.affineMap(norm_mat, []); % normalized state
        % step 2: compute brake output of rl_controller
        brake = rl_controller.reach(norm_X);
        % step 3: compose brake output and speed
        speed_brakes = get_speed_brakes(brake, norm_X);
        % step 4: compute transformer output
        tf_outs = get_tf_outs(transformer, speed_brakes);
        % step 5: get control output
        controls = get_controls(tf_outs, norm_X);
        % step 6: scale control output
        scaled_controls = scale_controls(controls, scale_mat);
        % step 7: compute reachable set of the plant with the new control input and
        % initial state
        X0 = stepReachPlant(A, B, X0, scaled_controls);  
        S = [S X0];
    end
end


function inv_TTC_max = get_inv_TTC_max(S)
    n = length(S); % number of stars in the reachable set
    B = [];
    for i=1:n
        if isa(S(i), 'Star')
            B = [B S(i).getBox]; % get box bounds the star set
        else
            B = [B S(i)];
        end
    end

    % Estimate minimum value of the function g = v^2 + 2ad
    g_min =  zeros(n,1);
    g_max = zeros(n,1);

    for i=1:n

        fun_min = @(x)x(2)^2 + 2*x(1)*x(3);
        fun_max = @(x)-x(2)^2 - 2*x(1)*x(3);
        lb = B(i).lb;
        ub = B(i).ub;
        x0 = lb;   
        [~,fval] = fmincon(fun_min,x0,[],[],[],[],lb,ub);
        g_min(i) = fval;
        [~,fval] = fmincon(fun_max,x0,[],[],[],[],lb,ub);
        g_max(i) = -fval;
    end

    % compute TTC^-1
    inv_TTC_min = zeros(n,1);
    inv_TTC_max = zeros(n,1);

    for i=1:n

        if g_max(i) <= 0

            inv_TTC_min(i) = 0;
            inv_TTC_max(i) = 0;

        elseif g_min(i) >= 0

            fun_min = @(x)(-x(3)/(x(2) - sqrt(x(2)^2 + 2*x(1)*x(3))));
            fun_max = @(x)(x(3)/(x(2) - sqrt(x(2)^2 + 2*x(1)*x(3))));
            lb = B(i).lb;
            ub = B(i).ub;
            x0 = lb; 
            if lb(3) == 0 && ub(3) == 0
                inv_TTC_min(i) = 0;
                inv_TTC_max(i) = 0;
            else
                [~,fval] = fmincon(fun_min,x0,[],[],[],[],lb,ub);
                inv_TTC_min(i) = fval;
                [~,fval] = fmincon(fun_max,x0,[],[],[],[],lb,ub);
                inv_TTC_max(i) = -fval;
            end

        elseif g_min(i) < 0 && g_max(i) > 0

            nonlcon = @feasiblecon;
            fun_max = @(x)(x(3)/(x(2) - sqrt(x(2)^2 + 2*x(1)*x(3))));
            lb = B(i).lb;
            ub = B(i).ub;
            x0 = lb;
            j = 0;
            while (j < 1000)
                for k=1:3
                        x0(k) = (ub(k) - lb(k)).*rand(1, 1) + lb(k);
                end
                x1 = fun_max(x0);
                if isreal(x1)
                    break;
                end
            end
            [~,fval] = fmincon(fun_max,x0,[],[],[],[],lb,ub,nonlcon);
            inv_TTC_max(i) = -fval;
            inv_TTC_min(i) = 0;       

        end

    end
    
    inv_TTC_max = max(inv_TTC_max);

end

function [c, ceq] = feasiblecon(x)
    c = -x(2)^2 - 2*x(1)*x(3);
    ceq = [];
end


function speed_brakes = get_speed_brakes(brake, norm_X)
    n = length(brake);
    speed_brakes = [];
    
    for i=1:n
        V = [norm_X.V(2, :);brake(i).V];
        speed_brakes = [speed_brakes Star(V, brake(i).C, brake(i).d)];
        
    end 
end

function tf_outs = get_tf_outs(transformer, speed_brakes)
    n = length(speed_brakes);
    tf_outs = [];
    parfor i=1:n
        tf_out = transformer.reach(speed_brakes(i));
        tf_outs = [tf_outs tf_out];
    end
end

function controls = get_controls(tf_outs, norm_X)
    n = length(tf_outs);
    controls = [];
    for i=1:n
        V = [norm_X.V(2, :);tf_outs(i).V];
        controls = [controls Star(V, tf_outs(i).C, tf_outs(i).d)];        
    end
end

function scaled_controls = scale_controls(controls, scale_mat)
    n = length(controls);
    scaled_controls = [];
    for i=1:n
        scaled_controls = [scaled_controls controls(i).affineMap(scale_mat, [])];
    end
end

function X = stepReachPlant(A, B, X0, scaled_controls)
    
    n = length(scaled_controls);
    X = X0.affineMap(A, []);
    U = [];
    for i=1:n
        U = [U scaled_controls(i).affineMap(B, [])];
    end
    
    next_X = [];
    for i=1:n
        V = X.V + U(i).V;
        next_X = [next_X Star(V, U(i).C, U(i).d)];
    end
    
    X = Star.get_hypercube_hull(next_X);
    X = X.toStar;
    
end