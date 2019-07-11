classdef IPDPolyhedron < Polyhedron
    methods
        function obj = IPDPolyhedron(data)
            % Polyhedron implicitly defined by primal/dual feasibility conditions

            global MPTOPTIONS
            % CR = { x | A*primal<=b+pB*x, dual>=0 }
            %    = { x | A*(F*x+g)<=b+pB*x, M*x+m>=0 }
            %    = { x | (A*F-pB)*x<=b-A*g, -M*x<=m }
            %    = { x | H*x<=h }
            A = [ data.OptProb.A*data.Primal.F-data.OptProb.pB; -data.DualIneq.F ];
            b = [ data.OptProb.b-data.OptProb.A*data.Primal.g; data.DualIneq.g ];
            b = b + MPTOPTIONS.modules.solvers.enum_pqp.ineq_backoff; % backoff tuned w.r.t. test_enum_pqp_10
            % include Ath*x<=bth bounds
            A = [A; data.OptProb.Ath];
            b = [b; data.OptProb.bth];
            if data.OptProb.me>0
                Ae = data.OptProb.Ae*data.Primal.F-data.OptProb.pE;
                be = data.OptProb.be-data.OptProb.Ae*data.Primal.g;
                args = {'H', [A, b], 'He', [Ae, be]};
            else
                % cheaper construction via Polyhedron(A, b)
                args = {A, b};
            end
            obj = obj@Polyhedron(args{:}); % must be a top-level command!
            obj.Data = data;
        end
        
        function P = toPolyhedron(obj)
            % converts the set to Polyhedron object
            P = [];
            for i = 1:numel(obj)
                p = Polyhedron('H', obj(i).H, 'He', obj(i).He);
                p.copyFunctionsFrom(obj(i));
                if numel(P)==0
                    P = p;
                else
                    P(i) = p;
                end
            end
        end
        
        function display(obj)
            % Display for IPDPolyhedron
            
            if numel(obj)>1
                fprintf('Array of %d implicit polyhedra.\n', numel(obj));
            else
                fprintf('Implicit polyhedron in %dD.\n', obj.Dim);
            end
        end
        
        function tf = contains(P, x, varargin)
            % Tests x \in P
            
            global MPTOPTIONS
            if isempty(MPTOPTIONS)
                MPTOPTIONS = mptopt;
            end
            tolerance = MPTOPTIONS.abs_tol;
            
            narginchk(2, 3);
            if ~isnumeric(x)
                % containment of sets
                tf = contains@Polyhedron(P, x, varargin{:});
                return
            end
            
            tf = false(size(P));
            for i = 1:numel(P)
                dual_f = P(i).Data.DualIneq;
                dual = dual_f.F*x+dual_f.g;
                if isempty(dual), dual=zeros(1, P(i).Dim); end
                
                % check dual feasibility first
                if min(dual)>=-tolerance
                    primal_f = P(i).Data.Primal;
                    primal = primal_f.F*x+primal_f.g;
                    prob = P(i).Data.OptProb;
                    % check primal feasibility
                    if max(prob.A*primal-prob.b-prob.pB*x)<=tolerance
                        if prob.me>0
                            % also check equalities
                            if norm(prob.Ae*primal-prob.be-prob.pE*x)<=tolerance
                                tf(i) = true;
                            end
                        else
                            tf(i) = true;
                        end
                    end
                end
            end
        end
    end
       
end
