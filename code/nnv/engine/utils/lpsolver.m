function [fval, exitflag] = lpsolver(f, A, b, Aeq, Beq, lb, ub, lp_solver, opts)
    % common function for al linear programming problems
    % for now we will set linprog with default options to solve every task,
    % seems to be faster than glpk or yalmip, and if it fails, we try glpk
    % (glpk seems a little more robust than linprog, fails less)
    % But if users prefer other solvers, we can allow them to use those.
    % The other default one is glpk (installed with NNV setup), the other
    % option will be to use yalmip and let users specify the solver they
    % want to use and call yalmip from this function

    if ~exist('lp_solver', 'var') || ~isempty(Aeq) || ~isempty(Beq)
        lp_solver = 'linprog'; % default solver, sometimes fails with mpt, so glpk as backup
    end

    if ~exist('opts', 'var')
        opts = 'normal';
    elseif ~strcmp(opts, 'emptySet')
        opts = 'normal';
    end

    % ensure variables are all of type double
    f = double(f); A = double(A); b = double(b); lb = double(lb);
    Aeq = double(Aeq); Beq = double(Beq); ub = double(ub); 
    
    % Solve using linprog (glpk as backup)
    if strcmp(lp_solver, 'linprog')
        options = optimoptions(@linprog, 'Display','none'); 
        options.OptimalityTolerance = 1e-10; % set tolerance
        % first try solving using linprog
        [~, fval, exitflag, ~] = linprog(f, A, b, Aeq, Beq, lb, ub, options);
        % found (1), not feasible point found (-2), infeasible (-5)
        if ~( (exitflag == -2 && strcmp(opts, 'emptySet')) || ismember(exitflag, [1, -5]) ) % return existflag if one of these conditions is met
            if ~isempty(Aeq) || ~isempty(Beq)
                error("Problem cannot be solved by linprog, and task not supported by glpk.")
            else
                [~, fval, exitflag, ~] = glpk(f, A, b, lb, ub);
                if ~ismember(exitflag, [2, 5, 3, 4, 110]) % feasible (2), optimal (5), not feasible (3, 4, 110)
                    error("LP solver error. Task failed to be solved by linprog and glpk. GLPK exitflag = " + string(exitflag));
                end
                exitflag = "g" + string(exitflag); % e.g. g2 or g5
            end
        else
            exitflag = "l" + string(exitflag); % e.g. l1 
        end

    % solve using glpk (linprog as backup)
    elseif strcmp(lp_solver, 'glpk')
        [~, fval, exitflag, ~] = glpk(f, A, b, lb, ub);
        if ~ismember(exitflag, [2, 5, 3, 4, 110]) % feasible (2), optimal (5), not feasible (3, 4, 110)
            options = optimoptions(@linprog, 'Display','none'); 
            options.OptimalityTolerance = 1e-10; % set tolerance
            % first try solving using linprog
            [~, fval, exitflag, ~] = linprog(f, A, b, Aeq, Beq, lb, ub, options);
            if ~ismember(exitflag, [1,-2, -5]) % found (1), not feasible point found (-2), infeasible (-5)
                error("Problem cannot be solved by linprog, and task not supported by glpk. LINPROG exitflag = " + string(exitflag));
            end
            exitflag = "l" + string(exitflag); % e.g. l1 
        else
            exitflag = "g" + string(exitflag); % e.g. g5
        end

    else % Try using yalmip with defined solver, much slower in our tests
        ops = sdpsettings('solver',lp_solver, 'verbose', 0);
        x = sdpvar(length(stars(rs).predicate_lb),1);
        if isempty(Aeq) && isempty(Beq)
            if isempty(lb) && isempty(ub)
                constraints = A*x <= b;
            else
                constraints = [A*x <= b, lb <= x , x <= ub];
            end
        else
            constraints = [A*x <= b, Aeq*x == Beq, lb <= x , x <= ub];
        end
        diagnostics = optimize(constraints, f*x, ops);
        if diagnostics.problem == 0
            x = value(x);
            fval = f*x;
            exitflag = "l1"; % treat it as an optimal solution found with linprog
        else
            error('YALMIP solver thinks it is infeasible or it ran into an error.')
        end
        
    end
end