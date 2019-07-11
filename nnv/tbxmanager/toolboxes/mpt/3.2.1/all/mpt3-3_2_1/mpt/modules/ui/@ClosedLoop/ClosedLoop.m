classdef ClosedLoop < MPTUIHandle & IterableBehavior
    % Object representation of a closed-loop system
    %
    % Constructor:
    %   loop = ClosedLoop(controller, system)
    
    properties(SetAccess=private)
        controller % Controller
        system % Controlled system
    end
    
    methods
		
        function obj = ClosedLoop(controller, system)
            % Constructor:
            %   loop = ClosedLoop(controller, system)
            
            if nargin == 0
                return
			end
			
			% TODO: check compatibility of number of states/inputs
            assert(isa(system, 'AbstractSystem'), 'Invalid type of the second argument.');
            assert(controller.nu==system.nu, 'Incompatible number of inputs.');
            assert(controller.nx==system.nx, 'Incompatible number of states.');
            obj.system = system;
            obj.controller = controller;
		end
		
		function I = invariantSet(obj, varargin)
			% Computes invariant subset of the closed-loop system
			
			if ~obj.controller.isExplicit()
				error('Only explicit controllers supported.');
			end
			
			I = obj.toSystem().invariantSet(varargin{:});
		end
		
		function out = simulate(obj, x0, N_sim, varargin)
			%
			% Simulates the closed-loop system for a given number of steps
			%
			
			narginchk(3, Inf);
			error(validate_vector(x0, obj.system.nx, 'initial state'));

			% internal helper to derimine if "string" ends up with "part"
			function out = endswith(string, part)
				out = length(string)>=length(part) && ...
					isequal(string(end-length(part)+1:end), part);
			end

			% determine whether we have free references. if so, allow the
			% user to specify their respective full closed-loop profiles
			references = struct('value', {}, 'position', {});
			if nargin>3
				if mod(length(varargin), 2)~=0
					error('Options must come in key/value pairs.');
				end
				% find position of reference options in varargin
				ref_positions = find(cellfun(@(z) endswith(z, '.reference'), varargin));
				for i = 1:length(ref_positions)
					% references found, store their values and position of
					% the value in varargin such that we can updated them
					% later
					ref_name = varargin{ref_positions(i)};
					ref_value = varargin{ref_positions(i)+1};
					% validate dimensions: number of columns of ref_value
					% must either be 1 or Nsim 
					if size(ref_value, 2)==1
						% expand throughout N_sim
						ref_value = repmat(ref_value, 1, N_sim);
					elseif size(ref_value, 2)~=N_sim
						error('"%s" must have either 1 or %d columns.', ...
							ref_name, N_sim);
					end
					% store the reference
					ref.value = ref_value;
					% position of the reference value in varargin
					ref.position = ref_positions(i)+1;
					references(end+1) = ref;
				end
            end
			
            include_disturbance = isa(obj.system, 'ULTISystem');
            
			X = x0(:); U = []; Y = []; D = []; J = [];
			for k = 1:N_sim

				if k>1
					% update u.previous and y.previous
					u_position = find(cellfun(@(z) isequal(z, 'u.previous'), varargin));
					if ~isempty(u_position)
						varargin{u_position+1} = u;
					end
					y_position = find(cellfun(@(z) isequal(z, 'y.previous'), varargin));
					if ~isempty(y_position)
						varargin{y_position+1} = y;
					end
					% TODO: reject updating variables which are not
					% simulated (e.g. "d", "z")
				end
				
				% use the k-th step references
				for i = 1:length(references)
					varargin{references(i).position} = references(i).value(:, k);
				end
				
				[u, feasible, openloop] = obj.controller.evaluate(x0, varargin{:});
				if ~feasible
					warning('No control action found at step %d, aborting.', k);
					break
                end
                
                % note: we could use obj.system.update(u) here, but that
                % introduces significant overhead. update_equation is much
                % faster.
                if include_disturbance
                    [x0, y, ~, d] = obj.system.update_equation(x0, u);
                else
                    [x0, y] = obj.system.update_equation(x0, u);
                end
				X = [X x0];
				U = [U u];
				Y = [Y y];
                if include_disturbance
                    D = [D d];
                end
				J = [J openloop.cost];
			end
			out.X = X;
			out.U = U;
			out.Y = Y;
            if include_disturbance
                out.D = D;
            end
			out.cost = J;
		end

		function out = toSystem(obj)
			% Converts a closed loop system into a LTI/PWA system
			
			if ~obj.controller.isExplicit()
				error('Only explicit controllers supported.');
			end
			
			if numel(obj.controller.optimizer)>1
				error('Overlapping partitions not supported.');
			end
			
			if obj.controller.optimizer.Dim ~= obj.controller.model.nx
				error('Tracking controllers not supported.');
            end

            if isa(obj.system, 'ULTISystem') && length(obj.controller.optimizer.Set)==1
                % Uncertain LTI system + linear controller
                
                feedback = obj.controller.optimizer.Set.getFunction('primal');
                F = feedback.F(1:obj.system.nu, :);
                g = feedback.g(1:obj.system.nu);
                assert(nnz(g)==0, 'Only linear controllers are supported.');
                unc = obj.system.cellifyMatrices();

                % 1) domain of the closed-loop system
                A = []; b = [];
                % umin <= F*x+g <= umax + set constraints
                Bu = obj.system.u.boundsToPolyhedron();
                A = [A; Bu.A*F];
                b = [b; Bu.b-Bu.A*g];
                % xmin <= x <= xmax + set constraints
                Bx = obj.system.x.boundsToPolyhedron();
                A = [A; Bx.A];
                b = [b; Bx.b];
                % ymin <= (C*x+D*u) <= ymax
                By = obj.system.y.boundsToPolyhedron();
                for ic = 1:numel(unc.C)
                    for id = 1:numel(unc.D)
                        A = [A; By.A*(unc.C{ic}+unc.D{id}*F)];
                        b = [b; By.b-By.A*unc.D{id}*g];
                    end
                end
                D = Polyhedron('A', A, 'b', sanitize_inf(b));
                D = D.intersect(obj.system.domainx);

                % 2) uncertain dynamics of the closed-loop system
                nA = numel(unc.A);
                nB = numel(unc.B);
                An = cell(1, nA*nB);
                i = 0;
                for ia = 1:nA
                    for ib = 1:nB
                        i = i + 1;
                        An{i} = unc.A{ia}+unc.B{ib}*F;
                    end
                end
                out = ULTISystem('A', An, 'C', obj.system.C, ...
                    'D', obj.system.D, 'domain', D);
                out.d = obj.system.d.copy();
                
            elseif isa(obj.system, 'ULTISystem')
                error('Uncertain LTI systems with PWA controllers not yet supported.');
                
            elseif isa(obj.system, 'LTISystem') && length(obj.controller.optimizer.Set)==1
				% LTI system + linear controller = LTI system with
				% domain restricted to the set where the controller
				% satisfies constraints (not necessarily recursively,
				% though)
				
				feedback = obj.controller.optimizer.Set.getFunction('primal');
				F = feedback.F(1:obj.system.nu, :);
				g = feedback.g(1:obj.system.nu);
				
				% 1) domain of the closed-loop system
				A = []; b = [];
				
                % umin <= F*x+g <= umax + set constraints
                Bu = obj.system.u.boundsToPolyhedron();
                A = [A; Bu.A*F];
                b = [b; Bu.b-Bu.A*g];
                
                % xmin <= x <= xmax + set constraints
                Bx = obj.system.x.boundsToPolyhedron();
                A = [A; Bx.A];
                b = [b; Bx.b];
                
				% ymin <= y <= ymax + set constraints
                By = obj.system.y.boundsToPolyhedron();
                A = [A; By.A*(obj.system.C+obj.system.D*F)];
                b = [b; By.b-By.A*obj.system.D*g];
				
				D = Polyhedron('A', A, 'b', sanitize_inf(b));
				D = D.intersect(obj.system.domainx);
				
				% 2) dynamics of the closed-loop system
				An = obj.system.A + obj.system.B*F;
				Bn = zeros(obj.system.nx, 0);
				Cn = obj.system.C + obj.system.D*F;
				Dn = zeros(obj.system.ny, 0);
				fn = obj.system.B*g;
				gn = obj.system.D*g;
				
				% 3) construct the LTI system
				out = LTISystem('A', An, 'B', Bn, 'C', Cn, 'D', Dn, ...
					'f', fn, 'g', gn, 'domain', D);
				
			elseif isa(obj.system, 'LTISystem') || isa(obj.system, 'PWASystem')
				% LTI or PWA system + PWA controller = PWA system
				
				R = obj.controller.feedback.Set;
				Dom = obj.system.domain;
				
				A = {}; B = {}; C = {}; D = {}; f = {}; g = {};
				Rn = [];
				for i = 1:length(R)
					for j = 1:length(Dom)
                        % extract parameters of the affine control law u=F*x+G
                        F = R(i).Func{1}.F(1:obj.system.nu, :);
                        G = R(i).Func{1}.g(1:obj.system.nu);
                        % check intersection with the j-th dynamics
                        %   D_j = { [x; u] | H*[x; u] <= h }
                        Hx = Dom(j).A(:, 1:obj.system.nx);
                        Hu = Dom(j).A(:, obj.system.nx+1:end);
                        Dom_j = Polyhedron(Hx+Hu*F, Dom(j).b-Hu*G);
						P = R(i).intersect(Dom_j);
						if P.isFullDim
                            if iscell(obj.system.A)
                                sys_A = obj.system.A{j};
                                sys_B = obj.system.B{j};
                                sys_f = obj.system.f{j};
                                sys_C = obj.system.C{j};
                                sys_D = obj.system.D{j};
                                sys_g = obj.system.g{j};
                            else
                                sys_A = obj.system.A;
                                sys_B = obj.system.B;
                                sys_f = obj.system.f;
                                sys_C = obj.system.C;
                                sys_D = obj.system.D;
                                sys_g = obj.system.g;
                            end
                            A{end+1} = sys_A + sys_B*F;
							B{end+1} = zeros(obj.system.nx, 0);
							f{end+1} = sys_B*G + sys_f;
							C{end+1} = sys_C + sys_D*F;
							D{end+1} = zeros(obj.system.nu, 0);
							g{end+1} = sys_D*G + sys_g;
							Rn = [Rn, P];
						end
					end
				end
				out = PWASystem('A', A, 'B', B, 'C', C, 'D', D, ...
					'f', f, 'g', g, 'domain', Rn, 'Ts', obj.system.Ts);
                out.Internal.Domain = obj.controller.optimizer.Domain;
			else
				error('Unsupported system.');
			end
			out.x = obj.system.x.copy();
			out.u = SystemSignal;
			out.y = obj.system.y.copy();
		end
    end
end
