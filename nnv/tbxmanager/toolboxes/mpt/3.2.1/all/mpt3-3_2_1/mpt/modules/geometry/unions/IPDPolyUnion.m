classdef IPDPolyUnion < PolyUnion
    % Union of IPDPolyhedron sets
    methods
        function obj = IPDPolyUnion(sets)
            % Constructs union of IPDPolyhedron sets

            if nargin==0
                return
            end
            assert(numel(sets)>0, 'At least one set please.');
            for i = 1:numel(sets)
                % TODO: allow to import from Polyhedron/PolyUnion (probably
                % need to pass the optimization problem as an additional
                % input)
                assert(isa(sets(i), 'IPDPolyhedron'), 'The input must be an array of IPDPolyhedron objects.');
            end
            obj.Set = sets;
            obj.Dim = sets(1).Data.OptProb.d;
            obj.Data.OptProb = sets(1).Data.OptProb;
        end
        
        function toC(obj, filename, function_to_export, varargin)
            % Exports evaluation of a given function to a C code
            %
            %   obj.toC(filename, function_to_export)
            %
            % Once exported, the file "filename_mex.mex" will be
            % automatically generated. It can then be called by
            %   [Z, REGION] = filename_mex(X)
            % where X is the evaluation point, Z is the computed function
            % value, and REGION is the index of the "active" region. If X
            % is not in the domain of the function, Z=NaN and REGION=0 will
            % be returned.
            %
            % The function value Z is obtained by sequentially checking all
            % optimal active sets.
            %
            % The current implementation is pretty rudimetary. The C code
            % (in fact, the compiled mex) is obtained by running codegen
            % on the Matlab export obtained via obj.toMatlab().

            narginchk(3, Inf);
            % TODO: tailored C code generation without codegen
            mfname = obj.toMatlab(filename, function_to_export, varargin{:});
            fname = mfname(1:end-2); % strip the '.m' extension
            fprintf('Calling codegen on "%s"...\n', mfname);
            cmd = sprintf('codegen %s -args {zeros(%d,1)}', mfname, obj.Dim);
            eval(cmd);
            fprintf('Function "%s" exported to "%s_mex.%s"\n', function_to_export, fname, mexext);
        end
        
        function f_name = toMatlab(obj, filename, function_to_export, varargin)
            % Exports evaluation of a given function to a matlab code
            %
            %   obj.toMatlab(filename, function_to_export)
            %
            % Once exported, the file "filename.m" will become available.
            % It can then be called by
            %   [Z, REGION] = filename(X)
            % where X is the evaluation point, Z is the computed function
            % value, and REGION is the index of the "active" region. If X
            % is not in the domain of the function, Z=NaN and REGION=0 will
            % be returned.
            %
            % The function value Z is obtained by sequentially checking all
            % optimal active sets.
            %
            % If obj.toMatlab(fname, fun2export, 'directPrimal', false) is
            % used, which is the default, then primal variables are
            % computed from the duals. This leads to less data being
            % stored, but more on-line computations performed. If
            % directPrimal=true, we store the parametric primal optimizer.
            
            global MPTOPTIONS
            
            narginchk(2, Inf);
            p = inputParser;
            p.addOptional('tiebreak', 'first-region');
            p.addParamValue('tolerance', MPTOPTIONS.abs_tol);
            p.addParamValue('directPrimal', false);
            p.parse(varargin{:});
            options = p.Results;
            
            % TODO: support multiple unions
            assert(numel(obj)==1, 'Single union please.');
            assert(obj.hasFunction(function_to_export), 'No such function "%s" in the union.', function_to_export);
            
            % TODO: support generic tiebreaks
            assert(isequal(options.tiebreak, 'first-region'), 'Only the "first-region" tiebreak is supported for now.');
            % TODO: support problems with equalities
            assert(obj.Data.OptProb.me==0, 'Problems with equality constraints not yet supported.');

            % TODO: also export and check parametric bounds (Ath*x<=bth)
            
            % get canonical file name
            [~, f_name, ~] = fileparts(filename);
            % replace invalid characters by an underscore
            fun_name = regexprep(f_name, '[^a-zA-Z0-9_]', '_');
            f_name = [fun_name '.m'];
            fid = fopen(f_name, 'w');
            assert(fid>0, 'Couldn''t open file "%s" for writing.', f_name);
            c = onCleanup(@() fclose(fid));
            
            fprintf(fid, 'function [z,i] = %s(x) %%#codegen\n', fun_name);
            fprintf(fid, '%% [zopt, region] = %s(x)\n', fun_name);
            fprintf(fid, '%%\n');
            fprintf(fid, '%% codegen %s -args {zeros(%d,1)}\n', f_name, obj.Dim);
            fprintf(fid, 'assert(isreal(x));\n');
            fprintf(fid, 'assert(size(x,1)==%d);\n', obj.Dim);
            fprintf(fid, 'assert(size(x,2)==1);\n');
            fprintf(fid, 'xh=[x;1];\n');
            
            % constraints G*z<=w+E*x+tol
            fprintf(fid, '%% primal constraints: G*z<=w+E*x+tol\n');
            fprintf(fid, 'G=%s;\n', mat2str(obj.Data.OptProb.A));
            fprintf(fid, 'Ew=%s;\n', mat2str([obj.Data.OptProb.pB obj.Data.OptProb.b]));
            fprintf(fid, 'Ewxt=Ew*xh+%g;\n', options.tolerance);
            
            % duals and number of duals
            D = []; nD = [];
            for i = 1:obj.Num
                dual_fun = obj.Set(i).Data.DualIneq;
                duals = [dual_fun.F, dual_fun.g];
                if isempty(duals)
                    duals = zeros(1, obj.Set(i).Dim+1);
                end
                nas = size(duals, 1);
                nD = [nD; nas];
                D = [D; duals];
            end
            nD = cumsum([1; nD]); % because we have a variable number of duals
            fprintf(fid, '%% duals d=D*[x;1] and number of duals\n');
            fprintf(fid, 'D=%s;\n', mat2str(D));
            fprintf(fid, 'nD=%s;\n', mat2str(nD));
            
            if options.directPrimal
                % primals for KKT conditions
                P = [];
                for i = 1:obj.Num
                    primal_fun = obj.Set(i).Data.Primal;
                    P = [P; primal_fun.F, primal_fun.g];
                end
                nP = length(obj.Set(1).Data.Primal.g);
                fprintf(fid, '%% primals p=P*[x;1]\n');
                fprintf(fid, 'P=%s;\n', mat2str(P));
            else
                % recover primal via duals:
                % From optimality condition: H*z + pF*x + f + Ga'*La =  0
                % it follows that z = -inv(H)*(pF*x + f + Ga'*La).
                % But we store directly the inverse of the Hessian
                H = inv(obj.Data.OptProb.H);
                pFf=[obj.Data.OptProb.pF obj.Data.OptProb.f];
                fprintf(fid, '%% inverted Hessian\n');
                fprintf(fid, 'Hi=%s;\n', mat2str(H));
                fprintf(fid, '%% pF*x+f\n');
                fprintf(fid, 'pFf=%s;\n', mat2str(pFf));
                % need to store active sets
                AS = [];
                for i = 1:obj.Num
                    as = obj.Set(i).Data.ActiveSet;
                    if isempty(as), as = 0; end
                    AS = [AS; as(:)];
                end
                fprintf(fid, '%% active sets (indexed via nD)\n');
                fprintf(fid, 'A=%s;\n', mat2str(AS)); % the active set
                fprintf(fid, '%% inv(H)*(pF*x+f)\n');
                fprintf(fid, 'HpFfx=Hi*(pFf*xh);\n');
            end
            
            % function to be evaluated
            Z = [];
            for i = 1:obj.Num
                fun = obj.Set(i).Functions(function_to_export);
                % TODO: support generic functions
                assert(isa(fun, 'AffFunction'), 'Only affine functions can be exported for now.');
                Z = [Z; fun.F, fun.g];
            end
            nZ = length(fun.g);
            fprintf(fid, '%% function to evaluate z=Z*[x;1]\n');
            fprintf(fid, 'Z=%s;\n', mat2str(Z));
            
            % sequential search
            fprintf(fid, 'for i=1:%d\n', obj.Num);
            fprintf(fid, '\t%% compute duals\n');
            fprintf(fid, '\td=D(nD(i):nD(i+1)-1,:)*xh;\n');
            fprintf(fid, '\tif all(d>=%g)\n', -options.tolerance); % if all duals non-negative
            fprintf(fid, '\t\t%% all duals non-negative, compute primals\n');
            if options.directPrimal
                fprintf(fid, '\t\tp=P((i-1)*%d+1:i*%d,:)*xh;\n', nP, nP); % evaluate primals
            else
                fprintf(fid, '\t\ta=A(nD(i):nD(i+1)-1);\n');
                fprintf(fid, '\t\tif a(1)==0\n');
                %fprintf(fid, '\t\t\tp=-Hi*(pFfx);\n');
                fprintf(fid, '\t\t\tp=-HpFfx;\n');
                fprintf(fid, '\t\telse\n');
                %fprintf(fid, '\t\t\tp=-Hi*(pFfx+G(a,:)''*d);\n');
                fprintf(fid, '\t\t\tp=-HpFfx-Hi*(G(a,:)''*d);\n');
                fprintf(fid, '\t\tend\n');
                % recover primal from dual
            end
            fprintf(fid, '\t\tif all(G*p<=Ewxt)\n'); % if all primals feasible
            fprintf(fid, '\t\t\t%% all primals feasible\n');
            fprintf(fid, '\t\t\tz=Z((i-1)*%d+1:i*%d,:)*xh;\n', nZ, nZ); % evaluate function + return
            fprintf(fid, '\t\t\treturn\n');
            fprintf(fid, '\t\tend\n'); % endif all primals feasible
            fprintf(fid, '\tend\n'); % endif all duals non-negative
            fprintf(fid, 'end\n'); % end for
            fprintf(fid, 'z=NaN(%d,1);i=0;\n', nZ); % infeasible
            fprintf(fid, 'end\n'); % end function
            
            % fclose(fid); % no need, done by oncleanup
            fprintf('Function "%s" exported to "%s"\n', function_to_export, f_name);
            if nargout==0
                clear f_name
            end
        end
    end
end
