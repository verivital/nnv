function sol = project(obj, y)
% <matDoc>
% <funcName>Project</funcName>
% <shortDesc>Compute the projection of the point y onto this polyhedron.<longDesc/>
%
% <syntax>
% <desc>Compute the projection of the point y onto this polyhedron.
%
% Solve the optimization problem:
%
%   J(y) = min ||y-x||_2 s.t. x in Set
%
% and return an optimizer x.
%
% </desc>
% <input name='x'>Vector of length dim(d)</input>
% <output name='sol'>Structure with fields
%   x       - Projected point or [] if emptyset
%   flag    - mptopt return type
%   flagTxt - mptopt return type text
%   dist    - Distance from x to the set, or [] if emptyset
% </output>
% </syntax>
%
% </matDoc>

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% reject arrays
error(obj.rejectArray());

dim = obj.Dim;
if dim<1
    error('Cannot project with empty polyhedra.');
end
validate_realmatrix(y);
if ~isequal(size(y, 1), obj.Dim)
	error('Input argument must have %d rows.', obj.Dim);
end

%% Project points onto the polyhedra
n_points = size(y, 2);
sol(1, n_points) = struct('x',[],'exitflag',[],'dist',[]); 

for j = 1:n_points
    
    % (x-y)'(x-y) = x'x - 2*x'y + y'y
    qp   = obj.optMat;
 
    % semidefinite QP
    cost = obj.buildCost(-2*y(:,j), 2*eye(obj.Dim));
    qp.f = cost.f; qp.H = cost.H;    
	qp.quickqp = true;
    opt  = mpt_solve(qp);
    
    % if not feasible, retry with setting bounds on all variables - helps
    % to recover feasibility 
    if opt.exitflag ~= MPTOPTIONS.OK
       qp.lb(qp.lb<-MPTOPTIONS.infbound) =-MPTOPTIONS.infbound;
       qp.ub(qp.ub>MPTOPTIONS.infbound) = MPTOPTIONS.infbound;
       opt  = mpt_solve(qp);
    end
    
    sol(1,j).exitflag    = opt.exitflag;
    sol(1,j).dist    = Inf; % infeasible proble = infinite distance by convention
    if sol(1,j).exitflag == MPTOPTIONS.OK,
        sol(1,j).x       = opt.xopt(1:obj.Dim);
        sol(1,j).dist = norm(sol(1,j).x - y(:,j));
    end
end



end

