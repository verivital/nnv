function R = mpt_call_linprog(S)
%
% header file to be inserted from XML source

global MPTOPTIONS
if isempty(MPTOPTIONS) && ~S.test
    MPTOPTIONS = mptopt;
end

if ~strcmpi(S.problem_type,'LP')
    error('mpt_call_linprog: LINPROG solver does not solve %s problems!',S.problem_type);
end

% overwrite default settings
if S.test
    try
        % MOSEK
        % I really hate doing this. Hey, MOSEK, stop messing with linprog!
        options = mskoptimset('linprog');
        options.Display = 'off';
    catch
        options=optimset(optimset('linprog'),'Display','off');
    end
else
    options=MPTOPTIONS.modules.solvers.linprog;
end

% direct call to linprog
try
    [R.xopt,R.obj,exitflag,OUTPUT,R.lambda]=linprog(S.f,S.A,S.b,S.Ae,S.be,S.lb,S.ub,S.x0,options);
catch
    % sometimes linprog fails with an infeasibility error
    exitflag = 0;
end

if exitflag>0 
    R.how = 'ok';
    if S.test
        R.exitflag = 1;
    else
        R.exitflag = MPTOPTIONS.OK;
    end
else
    R.how = 'infeasible';
    if S.test
        R.exitflag = 2;
    else
        R.exitflag = MPTOPTIONS.INFEASIBLE;
    end
end

%R.lambda=[lambdav.ineqlin; lambdav.eqlin];

end
