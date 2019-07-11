function solvers=mpt_detect_solvers
% MPT_DETECT_SOLVERS Searches for solvers and test the functionality
%
%  usage: mpt_detect_solvers - force to scan for solvers and test if they
%  work
%
% structure of the code:
%    1, search for available solvers
%    2, check functionality (only for non-parametric solvers)
%    3, sort solvers according to priority number
%    4, return results


%% 1) solvers to search for {'NAME',{'file1','file2',...},preference}
% m-files have extension ".m" and mex-files are without extension
% preference - according to this number the solvers will be ordered
% from the smallest number (fastest solver) higher.
% These numbers are used during the first time initialization of the
% toolbox, the user can change it afterwards in "mptopt"
solvers_list = {...
    {'LINPROG',{'linprog.m'},100},...
    {'QUADPROG',{'quadprog.m'},90},...
    {'GLPK',{'glpkcc'},80}...
    {'CDD',{'cddmex'},70},...
    {'CLP',{'mexclp'},130},...
    {'QPOASES',{'qpOASES'},110},...
    {'GUROBI',{'gurobi'},20},...
    {'NAG',{'e04nf','e04mf','mexnaglp','mexnagqp'},40},...
    {'CPLEX',{'cplexint','Cplex'},30},...
    {'MOSEK', {'mosekopt'}, 75}, ...
    {'LCP',{'lcp'},50},...
    {'QPC',{'qpas','qpip'},60},...
    {'SEDUMI',{'sedumi.m'},140},...
    {'QPSPLINE',{'QPspline.m'},150},...
    };

% parametric solvers (not tested, only checked if they are on the path)
parametric_solvers_list = {...
    {'PLCP',{'mpt_plcp.m'},10},...
    {'MPQP',{'mpt_mpqp.m'},20},...
    {'MPLP',{'mpt_mplp.m'},30},...
    {'ENUMPLCP',{'mpt_enum_plcp.m'},40},...
    {'ENUMPQP',{'mpt_enum_pqp.m'},50},...
    {'RLENUMPQP',{'mpt_enum_pqp.m'},60},...
    };

parametric_solvers_types.LP = {'MPLP','PLCP','ENUMPLCP'};
parametric_solvers_types.QP = {'MPQP','PLCP','ENUMPLCP','ENUMPQP','RLENUMPQP'};
parametric_solvers_types.LCP = {'PLCP','ENUMPLCP'};

% found solvers as strings
disp('MPT searches for solvers on the path ...');
disp(' ');
found_solvers = SubFindSolvers(solvers_list);
parametric_found_solvers = SubFindSolvers(parametric_solvers_list);

if size(found_solvers,1)<1
    error('mpt_detect_solvers: No solvers found, cannot proceed.');
end

%% 2) test solvers for correctness
% common setup for a trivial problem
H = 1;
f = 1;
A = [1;-1];
b = [1;1];
lb = 0;
test = true;

LPsolvers = zeros(size(found_solvers,1),1);
QPsolvers = LPsolvers;
MILPsolvers = LPsolvers;
MIQPsolvers = LPsolvers;
LCPsolvers = LPsolvers;

for i=1:size(found_solvers,1)
    % pick up the solver
    solver = found_solvers{i,1};
    
    % fore licensed software CPLEX or GUROBI, we need to test the license
    if strcmpi(solver,'cplex')
        fprintf('Checking CPLEX license ... ');
        if isempty(getenv('ILOG_LICENSE_FILE'))
            disp('Did not find an environment variable ILOG_LICENSE_FILE with CPLEX license.');
            %error('mpt_detect_solvers: CPLEX license not found.')
        else
            fprintf('ok\n');
        end
    end
    if strcmpi(solver,'gurobi')
        fprintf('Checking GUROBI license ... ');
        if isempty(getenv('GRB_LICENSE_FILE'))
            disp('Did not find an environment variable GRB_LICENSE_FILE with GUROBI license.');
            %error('mpt_detect_solvers: GUROBI license not found.')
        else
            fprintf('ok\n');
        end
%             % since GUROBI license is limited for short period of time, we activate
%             % it whenever GUROBI solver is initialized
%             if exist('grbvalidate','file')==2
%                 % "which" function does not find "grbvalidate" but we
%                 % know that it is in the same directory as "gurobi_mex"
%                 p = which('gurobi_mex');
%                 % cut off "gurobi_mex" name from the last filesep
%                 fp = strfind(p,filesep);
%                 p(fp(end)+1:end) = [];
%                 p = [p, 'grbvalidate'];
%                 % try to run validation script
%                 status = system(p);
%                 if status~=0
%                     disp('Executing of a validation script failed. Continuing anyway ...')
%                 end
%             else
%                 disp('mpt_detect_solvers: Failed to run "grbvalidate" to reactivate the license.');
%                 %error('mpt_detect_solvers: GUROBI activation file "grbvalidate" is missing. Did you install GUROBI correctly?');
%             end
        
    end
    
    % test LP solvers
    S = struct('solver',solver,'test',test,'f',f,'A',A,'b',b,'lb',lb);
    try
        R = mpt_solve(S);
        if R.exitflag==1
            LPsolvers(i) = true;
        end        
    end
    
    % test QP solvers  
    S = struct('solver',solver,'test',test,'H',H,'f',f,'A',A,'b',b,'lb',lb);
    try
        R = mpt_solve(S);
        if R.exitflag==1
            QPsolvers(i) = true;
        end
    end

    
    % test MILP solvers
    S = struct('solver',solver,'test',test,'f',f,'A',A,'b',b,'lb',lb,'vartype','B');
    try
        R = mpt_solve(S);
        if R.exitflag==1
            MILPsolvers(i) = true;
        end
    end
    
    % test MIQP solvers
    S = struct('solver',solver,'test',test,'H',H,'f',f,'A',A,'b',b,'lb',lb,'vartype','B');
    try
        R = mpt_solve(S);
        if R.exitflag==1
            MIQPsolvers(i) = true;
        end
    end
    
    % test LCP solvers
    S = struct('solver',solver,'test',test,'q',[f; b],'M',[H A'; -A zeros(2)]);    
    try 
        R = mpt_solve(S);
        if R.exitflag==1
            LCPsolvers(i) = true;
        end
    end
end

%% 3) sort solvers according to speed
LPorder = sortrows(found_solvers(logical(LPsolvers),:),2);
QPorder = sortrows(found_solvers(logical(QPsolvers),:),2);
MILPorder = sortrows(found_solvers(logical(MILPsolvers),:),2);
MIQPorder = sortrows(found_solvers(logical(MIQPsolvers),:),2);
LCPorder = sortrows(found_solvers(logical(LCPsolvers),:),2);

% assign solvers to MPTOPTIONS structure
solvers.LP = LPorder(:,1);
solvers.QP = QPorder(:,1);
solvers.MILP = MILPorder(:,1);
solvers.MIQP = MIQPorder(:,1);
solvers.LCP = LCPorder(:,1);

% if there is no LCP solver installed, PLCP solver will not work
if isempty(solvers.LCP)
    fprintf(['\nNo LCP solver detected.\nIt is strongly recommended to install LCP solver to have\n'...
          'numerically reliable results from the parametric LCP solver.\n'...
          'MPT will continue without the parametric LCP solver.\n']);
    % remove PLCP solver from the list
    ts=strcmp('PLCP',parametric_found_solvers(:,1));
    parametric_found_solvers(ts,:) = [];
end

% sort parametric solvers according to type of problem they solve
solvers.parametric.LP = parametric_found_solvers(ismember(parametric_found_solvers(:,1),parametric_solvers_types.LP),1);
solvers.parametric.QP = parametric_found_solvers(ismember(parametric_found_solvers(:,1),parametric_solvers_types.QP),1);
solvers.parametric.LCP = parametric_found_solvers(ismember(parametric_found_solvers(:,1),parametric_solvers_types.LCP),1);

% exit if no LP solver is working
if isempty(LPorder)
    error('mpt_detect_solvers: No working LP solver found, cannot proceed!');
end
% exit if no QP solver is working
if isempty(QPorder)
    error('mpt_detect_solvers: No working QP solver found, cannot proceed!');
end
% exit if no parametric solvers found
if isempty(solvers.parametric.QP) || isempty(solvers.parametric.LP)
    error('mpt_detect_solvers: No parametric solver found, cannot proceed!');
end

end

%% SUBFUNCTIONS
function found_solvers = SubFindSolvers(solvers_list)
% looks on the path for solvers given in "solvers_list" and shows them on
% the desktop

found_solvers = [];
for i=1:length(solvers_list)
    solver_name = solvers_list{i}{1};
    solver_file = solvers_list{i}{2};
    solver_pref = solvers_list{i}{3};
    solver_index = [];
    for j=1:length(solver_file)
        % instead of "exist" command, we rely on "which" because class
        % methods might not be visible
        if ~isempty(which(solver_file{j}))
            solver_index = [solver_index j];
        end
    end
    % prepare output to be put on display
    if ~isempty(solver_index)
        str = '';
        for j=1:length(solver_file)
            if j<length(solver_file)
                str = [solver_file{j},', ',str];
            else
                str = [str,solver_file{j}];
            end
        end
        % show found solvers on the screen
        fprintf(' %s %s %s \n',solver_name,...
            repmat('.', 1, 60-length(str)-length(solver_name)), str);
        found_solvers = [found_solvers; {solver_name} {solver_pref}];
    end
end
end
