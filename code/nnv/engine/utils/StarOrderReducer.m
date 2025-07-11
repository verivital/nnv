function reducedStar = StarOrderReducer(star, maxIter, tol)
% REDUCESTARALTERNATING compresses a Star set using alternating bilinear optimization
% Inputs:
%   star: struct with fields c, V, C, d, lb, ub
%   maxIter: maximum iterations per direction (e.g., 20)
%   tol: convergence tolerance for direction change (e.g., 1e-4)
%
% Output:
%   reducedStar: struct with c, V, lb, ub, empty C, d

% Unpack
c = star.c;
V = star.V;
C = star.C;
d = star.d;
lb = star.lb;
ub = star.ub;
m = size(V,1);
n = size(V,2);

% Initialize outputs
W = zeros(m,m);
l = zeros(m,1);
u = zeros(m,1);

options = optimoptions('linprog','Display','none');

% For each principal direction
V_proj = V;
c_proj = c;

for k = 1:m
    % Random initial w
    w = randn(m,1);
    w = w / norm(w);
    
    iter = 0;
    w_change = inf;
    
    while iter < maxIter && w_change > tol
        iter = iter + 1;
        w_old = w;
        
        % LPs to find min/max projection
        f_max = w' * V_proj;
        f_min = -f_max;
        
        [~, max_val] = linprog(f_max, C, d, [], [], lb, ub, options);
        [~, min_val] = linprog(f_min, C, d, [], [], lb, ub, options);
        min_val = -min_val;
        
        % Recover alpha maximizing and minimizing
        [alpha_max, ~] = linprog(f_max, C, d, [], [], lb, ub, options);
        [alpha_min, ~] = linprog(f_min, C, d, [], [], lb, ub, options);
        
        x_max = V_proj * alpha_max;
        x_min = V_proj * alpha_min;
        
        % Update w
        delta = x_max - x_min;
        w = delta / norm(delta);
        
        % Convergence check
        w_change = norm(w - w_old);
    end
    
    % Save direction
    W(:,k) = w;
    
    % Compute final intervals
    f_max = w' * V_proj;
    f_min = -f_max;
    [~, max_val] = linprog(f_max, C, d, [], [], lb, ub, options);
    [~, min_val] = linprog(f_min, C, d, [], [], lb, ub, options);
    l(k) = -min_val;
    u(k) = max_val;
    
    % Project remaining V onto orthogonal complement
    P = eye(m) - w * w';
    V_proj = P * V_proj;
end

% Output star
reducedStar.c = c;
reducedStar.V = W;
reducedStar.lb = l;
reducedStar.ub = u;
reducedStar.C = [];
reducedStar.d = [];
end

