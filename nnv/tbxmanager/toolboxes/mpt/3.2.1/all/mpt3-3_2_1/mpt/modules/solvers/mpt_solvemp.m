function R = mpt_solvemp(S)
%
% a common gateway function for parametric solvers
%

% input arguments is one object of OPT class
narginchk(1, 1);

if ~isa(S,'Opt')
    error('mpt_solvemp: Input argument must be an "Opt" class.');
end

% check parameters
if S.d<1
    error('mpt_solvemp: Problem does not contain parametric inputs.');
end

% call specific solver
switch upper(S.solver)
           
    case {'PLCP'}
        
        % call PLCP solver
        R = mpt_call_plcp(S);

    case {'ENUMPLCP'}
        
        % call ENUMPLCP solver
        R = mpt_call_enum_plcp(S);
        
    case {'MPQP'}
        
        % call MPQP solver
        R = mpt_call_mpqp(S);
        
    case {'MPLP'}
        
        % call MPLP solver
        R = mpt_call_mplp(S);
        
    case {'ENUMPQP'}
        % enumeration-based pQP solver
        R = mpt_enum_pqp(S, struct('regions', true));

    case {'RLENUMPQP'}
        % enumeration-based pQP solver
        R = mpt_enum_pqp(S, struct('regions', false));

    otherwise
        
        % unknown solver
        error('mpt_solvemp: Solver not recognized.');
end
        
    
    
end
