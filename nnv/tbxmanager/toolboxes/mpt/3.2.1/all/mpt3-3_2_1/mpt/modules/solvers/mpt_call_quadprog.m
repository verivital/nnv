function R = mpt_call_quadprog(S)
%
% header file to be inserted from XML source

global MPTOPTIONS
if isempty(MPTOPTIONS) && ~S.test
    MPTOPTIONS = mptopt;
end

if ~any(strcmpi(S.problem_type,{'QP'}))
    error('mpt_call_quadprog: QUADPROG solver does not solve %s problems!',S.problem_type);
end

% overwrite default settings
if S.test
    try
        % MOSEK
        % I really hate doing this. Hey, MOSEK, stop messing with linprog!
        options = mskoptimset('quadprog');
        options.Display = 'off';
    catch
        options=optimset(optimset('quadprog'),'Display','off');
    end
else
    options=MPTOPTIONS.modules.solvers.quadprog;
end

% make sure the Hessian is symmetric.options=optimset(optimset('quadprog'),'Display','off','LargeScale','off',...
if norm(S.H-S.H', Inf) < 1e-10,
    % but we only remove numerical noise if the hessian is only
    % "slightly" wrong. for clearly non-symmetrical hessians we still
    % let quadprog to display a proper warning
    S.H = (S.H + S.H')*0.5;
end

% direct call to quadprog
[R.xopt,R.obj,exitflag,OUTPUT,R.lambda]=quadprog(S.H,S.f,S.A,S.b,S.Ae,S.be,S.lb,S.ub,S.x0,options);
if exitflag>0 %then QUADPROG converged with a solution X.
    R.how = 'ok';
    if S.test
        R.exitflag = 1;
    else
        R.exitflag = MPTOPTIONS.OK;
    end
else
    % ==0 then the maximum number of iterations was exceeded (only occurs
    %     with large-scale method).
    % < 0 then the problem is unbounded, infeasible, or
    %     QUADPROG failed to converge with a solution X.
    R.how = 'infeasible';
    if S.test
        R.exitflag = 2;
    else
        R.exitflag = MPTOPTIONS.INFEASIBLE;
    end
end
%R.lambda=[lambdav.ineqlin; lambdav.eqlin];

end
