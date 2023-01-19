%% What is faster, glpk, linprog or yalmip?
% Note: Yalmip is not a solver, but a wrapper for a collection of solvers

% linprog syntax
% [~, fval, exitflag, ~] = linprog(f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub, options); 
% the empty args correspond to equalities (not really used for the star methods in NNV)

% glpk syntax
% [~, fval, exitflag, ~] = glpk(f, obj.C, obj.d, obj.predicate_lb, obj.predicate_ub); % no equalities here

% yalmip syntax (to solve the same problem)
% objective = f;
% x = sdpvar;
% contraints = [obj.C*x <= obj.d, lb <= x <= ub];
% diagnostics = optimize(constraints, objective) % possible 3rd argument to specify solver options
% fval = value(x);

s1 = load("TestStar1.mat");
s2 = load("TestStar3.mat");
s3 = load("TestStar4.mat");

stars = [s1.Rstar; s2.Rstar; s3.Rstar];

[time, fvals, flags] = test_solvers(stars);




%% Helper functions
function [time, fvals, flags] = test_solvers(stars)
    % Setup evaluation
    number_of_solvers = 3; % number of solvers
    dim = 5; % dimensionality for all stars
    number_of_stars = length(stars);
    % linprog specific
    linprogOptions = optimoptions(@linprog, 'Display','none');
    linprogOptions.OptimalityTolerance = 1e-10; % set tolerance
    % yalmip specific
    ops = sdpsettings('solver','linprog', 'verbose', 0);
    ops2 = sdpsettings('solver','mosek','verbose',1);
    ops3 = sdpsettings('solver','glpk', 'verbose', 1);

    % Memory allocation
    time = zeros(dim*number_of_stars, number_of_solvers); % timing
    flags = zeros(dim*number_of_stars, number_of_solvers); % exitflags
    fvals = zeros(dim*number_of_stars, number_of_solvers); % value computed

    % Start evaluating performance
    for rs = 1:number_of_stars
        for i=1:dim
            f = stars(rs).V(i, 2:stars(rs).nVar + 1); % objective function to solve for
            % linprog
            t = tic;
            [~, fval, exitflag, ~] = linprog(f, stars(rs).C, stars(rs).d, [], [], stars(rs).predicate_lb, stars(rs).predicate_ub, linprogOptions);
            t = toc(t);
            % results
            time((rs-1)*i+i,1) = t; 
            fvals((rs-1)*i+i,1) = fval; 
            flags((rs-1)*i+i,1) = exitflag;

            % glpk
            t = tic;
            [~, fval, exitflag, ~] = glpk(f, stars(rs).C, stars(rs).d, stars(rs).predicate_lb, stars(rs).predicate_ub);
            t = toc(t);
            % results
            time((rs-1)*i+i,2) = t; 
            fvals((rs-1)*i+i,2) = fval; 
            flags((rs-1)*i+i,2) = exitflag;
            
            % yalmip
            t = tic;
            x = sdpvar(length(stars(rs).predicate_lb),1);
            constraints = [stars(rs).C*x <= stars(rs).d, stars(rs).predicate_lb <= x , x <= stars(rs).predicate_ub];
            diagnostics = optimize(constraints, f*x, ops3);
            x = value(x);
            t = toc(t);
            fval = f*x;
            % results
            time((rs-1)*i+i,3) = t; 
            fvals((rs-1)*i+i,3) = fval; 
            flags((rs-1)*i+i,3) = diagnostics.problem;
        end
    end
end