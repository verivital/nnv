classdef Polyhedron < ConvexSet
	%<matDoc>
	% <funcName>Polyhedron</funcName>
	% <shortDesc>Create a polyhedron object.</shortDesc>
	% <longDesc>
	% Define a polyhedron:
	%
	%  P = {x | H*[x;-1] &lt;= 0} cap {x | He*[x;-1] = 0}
	%
	% or
	%
	%  P = {V'lam | lam &gt;= 0, sum(lam) = 1} + {R'gam | gam &gt;= 0}
	% </longDesc>
	%
	% <!-- The standard calling form -->
	% <syntax>
	%  <input  name='H' type='paramValue'>Inequality description (must be full-dimensional) {x | H*[x;-1] &lt;= 0}</input>
	%  <input  name='He' type='paramValue'>Affine set {x | He*[x;-1] = 0}</input>
	%  <input  name='V' type='paramValue'>Points defining the set conv(V)</input>
	%  <input  name='R' type='paramValue'>Rays defining the set cone(R)</input>
	%  <input  name='irredundantVRep' default='false' type='paramValue'>If true, then the given V-Rep is assumed irredundant</input>
	%  <input  name='irredundantHRep' default='false' type='paramValue'>If true, then the given H-Rep is assumed irredundant</input>
	%  <input  name='lb' type='paramValue'>Add a constraint of the form lb &lt; x</input>
	%  <input  name='ub' type='paramValue'>Add a constraint of the form x &lt; ub</input>
	%  <output name='P' class='Polyhedron'>The polyhedron</output>
	% </syntax>
	%
	% <!-- Copy constructor -->
	% <syntax>
	%  <desc>Copy constructor</desc>
	%  <input  name='Pin' class='Polyhedron'>Polyhedron to copy</input>
	%  <output name='P' class='Polyhedron'>Copy of Pin</output>
	% </syntax>
	%
	% <!-- Vertex form -->
	% <syntax>
	%  <desc>Shorthand for defining Vrep polyhedra</desc>
	%  <input  name='V'>Vertices of polyhedron (row-wise)</input>
	%  <output name='P' class='Polyhedron'>Vrep Polyedron conv(V)</output>
	% </syntax>
	%
	%
	% <!-- Inequality form -->
	% <syntax>
	%  <desc>Shorthand for defining Hrep polyhedra</desc>
	%  <input  name='A'>Normals of constraints (row-wise)</input>
	%  <input  name='b'>Offsets of constraints (row-wise) length(b)==size(A,1)</input>
	%  <output name='P' class='Polyhedron'>Hrep Polyhedron {x | Ax &lt;=b}</output>
	% </syntax>
	%
	% <syntax>
	%  <desc>Convert from YALMIP polyhedron</desc>
	%  <input  name='F'>Yalmip constraints</input>
	%  <input  name='vars'>All variables involved in F</input>
	%  <output name='P' class='Polyhedron'>Hrep Polyhedron {x | Ax &lt;=b}</output>
	% </syntax>
	%
	% <seeAlso>affineMap</seeAlso>
	% <seeAlso>contains</seeAlso>
	%
	% </matDoc>
	%
	
	
	% A polyhedron.
	%
	% This class repreents the following polyhedra:
	%  - unbounded polyhedra
	%  - lower-dimensional polyhedra
	%  - non-pointed polyhedra, or those that do not contain any vertices
	%  (e.g. <eq>\{x\in\mathbb{R}^2\,|\,x_1 \le 0\}</eq>)
	%
	% We examine each of these three points in the following sections.
	% <p>
	% <b>Unbounded Polyhedra</b><p>
	% All polyhedra can be decomposed into the sum of a bounded polytope a cone:
	% <eq>  P = \operatorname{conv}(V) + \operatorname{cone}(R)</eq><br>
	% and satisfy the Minkowski-Weyl theorem and can therefore be represented
	% either as the
	% intersection of a finite number of inequalities, or as the convex
	% combination of a finite
	% number of vertices (or rays).
	% <p>
	% <i>MPT3.0 will store all irredundant polyhedra as a decomposition into a polytope and a cone.</i><p>
	% <p>
	% <b>Lower-dimensional Polyhedra</b><p>
	% Theoretically there is no difference between full-dimensional and lower-dimensional
	% polyhedra either in representation or in the algorithms that operate on them. However,
	% experience has shown that if the affine hull of the polyhedron is not taken into account
	% explicitly, then virtually all algorithms will fail.
	% <p>
	% <i>MPT3.0 will store a polyhedron as the intersection of a full-dimensional polyhedron and an affine hull.</i>
	% <p>
	% <b>Unpointed Polyhedra</b><p>
	% Operations on polyhedra with non-empty lineality space (i.e. unpointed polyhedra) adds
	% significant complexity and difficulty. Thankfully, all convex sets can be decomposed into
	% the Minkowski sum of their lineality space, with their restriction onto a linear subspace
	% U that is perpendicular to <eq>\operatorname{lineal}(P)</eq><br>
	% <eq>P= \operatorname{lineal}(P) + (P \cap U)</eq><br>
	% where <eq>P \cap U</eq> has an empty lineality space. Therefore, it is always possible to
	% represent an unpointed polyhedron as the lifting of a pointed one:<br>
	% <eq>P = \{Ey\,|\,y \in \tilde P\}</eq><br>
	% where <eq>\tilde P \in \mathbb{R}^m</eq> is pointed and <eq>E \in \mathbb{R}^{n\times m}</eq> with n > m.
	% <p>
	% The requirement of dealing with an additional lifting operation for all polyhedra will add
	% more coding complexity to \mptThree than desired. Therefore, in the interest of
	% simplicity, MPT3.0 will handle unpointed polyhedra only in halfspace form, where the
	% lifting map is not required, and any operation that cannot function on lifted polyhedra
	% will return an error.
	% <p>
	% <i>MPT3.0 will have a limited ability to operate on unpointed polyhedra.  It will
	% handle unpointed polyhedra only in halfspace form, where the lifting map is not
	% required, and any operation that does not function on unpointed
	% polyhedra will throw an exception.</i>
	%
	% @author cjones
	%
	% see also Polyhedron.Polyhedron
	
	
	properties(SetAccess=private, GetAccess=private, Hidden)
		H_int  = []; % Inequality description {x | H*[x;-1] <= 0}
		He_int = []; % Affine set description {x | He*[x;-1] = 0}
		V_int  = []; % Vertex description {V'*lam | lam >= 0, sum(lam) = 1}
		R_int  = []; % Cone description {R'*gam | gam >= 0}
	end
	
	properties(SetAccess = private)
		irredundantVRep = false; % True if the V-Representation is irredundant
		irredundantHRep = false; % True if the H-Representation is irredundant
		hasHRep = false; % True if the H-representation is available
		hasVRep = false; % True if the V-representation is available
	end
	
	properties(SetAccess=private, Transient=true, Hidden)
		% Pre-computed matrices to make solving optimization problems faster
		optMat = []
	end
	
	properties(Constant=true, Hidden)
		% version of the polyhedron object is used to do proper importing in
		% loadobj()
		%
		% History:
		%   1.0: September 25, 2012 (introduction of lazy getters)
		version = 1.0
	end
	
	properties(Dependent=true, SetAccess=private, Transient=true)
		A  % Inequality description { x | A*x <= b }
		b  % Inequality description { x | A*x <= b }
		Ae % Affine set description { x | Aeq*x == be }
		be % Affine set description { x | Aeq*x == be }
		H  % Inequality description { x | H*[x; -1] <= 0 }
		He % Affine set description { x | He*[x; -1] == 0 }
		R  % Rays of the polyhedron
		V  % Vertices of the polyhedron
	end
	
	methods
		% lazy getters
		function A = get.A(obj)
			if ~obj.hasHRep && obj.hasVRep
				obj.computeHRep();
			end
			A = obj.H_int(:, 1:end-1);
		end
		function b = get.b(obj)
			if ~obj.hasHRep && obj.hasVRep
				obj.computeHRep();
			end
			b = obj.H_int(:, end);
		end
		function Ae = get.Ae(obj)
			if ~obj.hasHRep && obj.hasVRep
				obj.computeHRep();
			end
			Ae = obj.He_int(:, 1:end-1);
		end
		function be = get.be(obj)
			if ~obj.hasHRep && obj.hasVRep
				obj.computeHRep();
			end
			be = obj.He_int(:, end);
		end
		function H = get.H(obj)
			if ~obj.hasHRep && obj.hasVRep
				obj.computeHRep();
			end
			H = obj.H_int;
		end
		function He = get.He(obj)
			if ~obj.hasHRep && obj.hasVRep
				obj.computeHRep();
			end
			He = obj.He_int;
		end
		function V = get.V(obj)
			if ~obj.hasVRep && obj.hasHRep
				obj.computeVRep();
			end
			V = obj.V_int;
		end
		function R = get.R(obj)
			if ~obj.hasVRep && obj.hasHRep
				obj.computeVRep();
			end
			R = obj.R_int;
		end
		function M = get.optMat(obj)
			M = obj.optMat;
			if isempty(M)
				M = obj.buildSetRepresentation;
				obj.optMat = M;
			end
		end
	end
	
	methods(Access = public)
		
		function obj = Polyhedron(varargin)
			% POLYHEDRON Create a polyhedron object.
			%
			% -------------------------------------------------------------------
			% Description
			% -------------------------------------------------------------------
			%
			% Define a polyhedron:
			%
			%  P = {x | H*[x;-1] <= 0} cap {x | He*[x;-1] = 0}
			%
			% or
			%
			%  P = {V'lam | lam >= 0, sum(lam) = 1} + {R'gam | gam >= 0}
			%
			% -------------------------------------------------------------------
			% Syntax
			% -------------------------------------------------------------------
			%
			% P = Polyhedron(varargin)
			%
			% Parameters specified as param/value pairs: (all optional)
			%  H  - Inequality description (must be full-dimensional) {x | H*[x;-1] <= 0}
			%  He - Affine set {x | He*[x;-1] = 0}
			%  V  - Points defining the set conv(V)
			%  R  - Rays defining the set cone(R)
			%  irredundantVRep - [false] If true, then the given V-Rep is
			%                    assumed irredundant
			%  irredundantHRep - [false] If true, then the given H-Rep is
			%                    assumed irredundant
			%  lb - Lower bound
			%  ub - Upper bound
			%  Data - arbitraty user data
			%
			% -------------------------------------------------------------------
			% Shorthand syntax
			%
			%  P = Polyhedron(Polyhedron). Copy constructor
			%  P = Polyhedron(V),    V = matrix. Create Vrep Polyhedron
			%  P = Polyhedron(A, b), A = matrix, b = vector. Create Hrep Polyhedron
			
			global MPTOPTIONS
			if isempty(MPTOPTIONS)
				MPTOPTIONS = mptopt;
			end
			
			obj.Functions = containers.Map;
			% Shorthand syntax:
			if nargin==0
				% empty Polyhedron object
				obj.Dim = 0;
				obj.H_int = zeros(0, 1);
				obj.He_int = zeros(0, 1);
				obj.Internal=struct('Empty',[],'Bounded',[],'FullDim',[],...
					'InternalPoint',[],'ChebyData',[]);
				return
				
			elseif nargin==2
				A = varargin{1}; b = varargin{2};
				if isnumeric(A) && isnumeric(b)
					% Hrep polyhedron without equality constraints: Polyhedron(A, b)
					[n, d] = size(A);
					if n~=length(b)
						error('Number of rows does not hold between arguments "A", "b".')
					end

					% replace a'*x<=+/-Inf by 0'*x<=+/-1
					InfRows = isinf(b);
					if any(InfRows)
						% normalize a'*x <= +/-Inf to 0'*x <= +/- 1
						A(InfRows, :) = 0*A(InfRows, :);
						b(InfRows) = sign(b(InfRows));
                    end

                    H = [A, b];
					% replace nearly-zero entries by zero
                    if issparse(H)
                        % much faster version of H(abs(H)<tol)=0 for sparse
                        % matrices (by J. Loefberg)
                        [i,j,k] = find(H);
                        k(abs(k)<MPTOPTIONS.zero_tol)=0;
                        H = sparse(i,j,k,size(H,1),size(H,2));
                    else
                        H(abs(H)<MPTOPTIONS.zero_tol) = 0;
                    end

					obj.Dim = d;
					obj.H_int = full(H);
					obj.He_int = zeros(0, d+1);
					obj.hasHRep = ~isempty(obj.H_int);
					obj.V_int = zeros(0, d);
					obj.R_int = zeros(0, d);
					obj.Internal=struct('Empty',[],'Bounded',[],'FullDim',[],...
						'InternalPoint',[],'ChebyData',[]);
					obj.Data = [];
					obj.irredundantVRep = false;
					obj.irredundantHRep = false;
					return
					
				elseif isa(A, 'lmi') && isa(b, 'sdpvar')
					% Convert from YALMIP format
					F = A;
					x = b;
					
					h = 0; % Cost func
					solving_parametric = 0;
					options = sdpsettings;
					[interfacedata,recoverdata,solver,diagnostic,F,Fremoved] = compileinterfacedata(F,[],[],h,options,0,solving_parametric);
					
					vars = find(ismember(recoverdata.used_variables,getvariables(x)));
					interfacedata.requested_variables = vars;
					mat = yalmip2mpt(interfacedata);
					
					% Often end up with inf in the RHS of YALMIP objects
					i = isinf(mat.W);
					mat.G(i,:) = [];
					mat.W(i,:) = [];
					obj = Polyhedron('H', [mat.G mat.W], 'He', [mat.Aeq mat.beq]);
					return
				end
				
			elseif nargin==1 && isa(varargin{1}, 'Polyhedron')
				% copy constructor
				obj = varargin{1}.copy();
				return
				
			elseif nargin==1 && isnumeric(varargin{1})
				% Vrep polyhedron: Polyhedron(V)
				obj = Polyhedron('V', varargin{1});
				return
				
			end
			
			% otherwise parse all inputs
			ip = inputParser;
			ip.KeepUnmatched = false;
			ip.addOptional('P', [], @validate_polyhedron);
			ip.addParamValue('H',  [], @validate_realmatrix);
			ip.addParamValue('He', [], @validate_realmatrix);
			ip.addParamValue('V',  [], @validate_realmatrix );
			ip.addParamValue('R',  [], @validate_realmatrix);
			ip.addParamValue('irredundantVRep', false, @validate_logicalscalar );
			ip.addParamValue('irredundantHRep', false, @validate_logicalscalar );
			ip.addParamValue('lb', [], @validate_realinfvector );
			ip.addParamValue('ub', [], @validate_realinfvector );
			ip.addOptional('Data', []); % user data can be anything
			ip.addOptional('A', [], @validate_realmatrix);
			ip.addOptional('b', [], @validate_realinfvector);
			ip.addOptional('Ae', [], @validate_realmatrix);
			ip.addOptional('be', [], @validate_realvector);
			
			ip.parse(varargin{:});
			p = ip.Results;
			
			% collect data from A, b, Ae, be if provided
			if ~isempty(p.A) && ~isempty(p.b)
				if size(p.A,1)~=numel(p.b(:))
					error('Number of rows does not hold between arguments "A", "b".')
				end
				if isempty(p.H)
					p.H =[p.A, p.b(:)];
				else
					if size(p.H,2)~=size(p.A,2)+1
						error('Number of columns does not hold between arguments "H", "A", "b".')
					end
					p.H = [p.H; p.A, p.b(:)];
				end
			end
			if ~isempty(p.Ae) && ~isempty(p.be)
				if size(p.Ae,1)~=numel(p.be(:))
					error('Number of rows does not hold between arguments "Ae", "be".')
				end
				if isempty(p.He)
					p.He = [p.Ae, p.be(:)];
				else
					if size(p.He,2)~=size(p.Ae,2)+1
						error('Number of columns does not hold between arguments "H", "A", "b".')
					end
					p.He = [p.He; p.Ae, p.be(:)];
				end
			end
			
			% Check that dimensions are correct
			d = max([size(p.H,2)-1 size(p.He,2)-1 size(p.V,2) size(p.R,2)]);
			%       if d == 0, error('Dimension must be strictly positive'); end
			if size(p.H,2)  > 0 && size(p.H,2)-1  ~= d || ...
					size(p.He,2) > 0 && size(p.He,2)-1 ~= d || ...
					size(p.V,2)  > 0 && size(p.V,2)    ~= d || ...
					size(p.R,2)  > 0 && size(p.R,2)    ~= d
				error('Input matrices must have the same dimension');
			end
			if d > 0
				if ~isempty(p.lb) && length(p.lb)   ~= d || ...
						~isempty(p.ub) && length(p.ub)  ~= d
					error('Upper lower bounds must be of length %i', d);
				end
			else
				d = max([length(p.lb) length(p.ub)]);
				if ~isempty(p.lb) && length(p.lb) ~= d || ...
						~isempty(p.ub) && length(p.ub) ~= d
					error('Upper lower bounds must be of length %i', d);
				end
			end
			obj.Dim = d;
			
			% check if every LB is actually lower than UB
			if ~isempty(p.lb) && ~isempty(p.ub)
				for i=1:d
					if p.lb(i)>p.ub(i) + MPTOPTIONS.zero_tol
						error('Polyhedron: Lower bound at element %d must not be greater than its upper bound.',i);
					end
				end
			end
			
			% Add upper/lower bounds
			if ~isempty(p.lb), p.H = [p.H;-eye(d) -p.lb(:)]; end
			if ~isempty(p.ub), p.H = [p.H; eye(d)  p.ub(:)]; end

			if ~isempty(p.H)
				% replace a'*x<=+/-Inf by 0'*x<=+/-1
				InfRows = isinf(p.H(:, end));
				if any(InfRows)
					% normalize a'*x <= +/-Inf to 0'*x <= +/- 1
					A = p.H(:, 1:end-1);
					b = p.H(:, end);
					A(InfRows, :) = 0*A(InfRows, :);
					b(InfRows) = sign(b(InfRows));
					p.H = [A, b];
				end
				
				% replace nearly-zero entries by zero
                if issparse(p.H)
                    % much faster version of H(abs(H)<tol)=0 for sparse
                    % matrices (by J. Loefberg)
                    [i,j,k] = find(p.H);
                    k(abs(k)<MPTOPTIONS.zero_tol)=0;
                    p.H = sparse(i,j,k,size(p.H,1),size(p.H,2));
                else
                    p.H(abs(p.H)<MPTOPTIONS.zero_tol) = 0;
                end
			end
			if ~isempty(p.He)
				% replace nearly-zero entries by zero
                if issparse(p.He)
                    % much faster version of H(abs(H)<tol)=0 for sparse
                    % matrices (by J. Loefberg)
                    [i,j,k] = find(p.He);
                    k(abs(k)<MPTOPTIONS.zero_tol)=0;
                    p.He = sparse(i,j,k,size(p.He,1),size(p.He,2));
                else
                    p.He(abs(p.He)<MPTOPTIONS.zero_tol) = 0;
                end
			end
			if ~isempty(p.V)
				% replace nearly-zero entries by zero
                if issparse(p.V)
                    % much faster version of H(abs(H)<tol)=0 for sparse
                    % matrices (by J. Loefberg)
                    [i,j,k] = find(p.V);
                    k(abs(k)<MPTOPTIONS.zero_tol)=0;
                    p.V = sparse(i,j,k,size(p.V,1),size(p.V,2));
                else
                    p.V(abs(p.V)<MPTOPTIONS.zero_tol) = 0;
                end
			end
			if ~isempty(p.R)
				% replace nearly-zero entries by zero
                if issparse(p.R)
                    % much faster version of H(abs(H)<tol)=0 for sparse
                    % matrices (by J. Loefberg)
                    [i,j,k] = find(p.R);
                    k(abs(k)<MPTOPTIONS.zero_tol)=0;
                    p.R = sparse(i,j,k,size(p.R,1),size(p.R,2));
                else
                    p.R(abs(p.R)<MPTOPTIONS.zero_tol) = 0;
                end
			end

			% Assign data
			obj.H_int  = full(p.H);
			obj.He_int = full(p.He);
			obj.V_int  = full(p.V);
			obj.R_int  = full(p.R);
			obj.irredundantVRep = false;
			obj.irredundantHRep = false;
			obj.hasHRep = ~isempty(p.H) || ~isempty(p.He);
			obj.hasVRep = ~isempty(p.V) || ~isempty(p.R);
			
			obj.Internal=struct('Empty',[],'Bounded',[],'FullDim',[],'InternalPoint',[],'ChebyData',[]);
			obj.Data = p.Data;
						
			if size(obj.R_int,1) > 0 || size(obj.V_int,1) > 0
				obj.irredundantVRep = p.irredundantVRep;
			end
			if size(obj.H_int,1) > 0
				obj.irredundantHRep = p.irredundantHRep;
			end
			
			% Fix sizes of empty data
			if isempty(obj.H_int),  obj.H_int  = zeros(0,d+1); end
			if isempty(obj.He_int), obj.He_int = zeros(0,d+1); end
			if isempty(obj.V_int),  obj.V_int  = zeros(0,d);   end
			if isempty(obj.R_int),  obj.R_int  = zeros(0,d);   end
			
			% Add the origin if this is a cone
			if size(obj.R_int,1)>0 && size(obj.V_int,1)==0
				% fprintf('Adding the origin to this cone');
				obj.V_int = zeros(1,d);
            end

			% Compute a minimum representation for the affine set
            if ~isempty(obj.He_int)
                if norm(obj.He_int)==0 && ...
                        (isempty(obj.H_int) || norm(obj.H_int)==0)
                    % corner case: 0*x=0
                    obj = Polyhedron.fullSpace(obj.Dim);
                else
                    % compute the minimal affine representation
                    obj.He_int = mpt_minAffineRep(obj.He_int);
                end
            end
			
		end

		function answer = isPointed(P)
			% returns true if the polyhedron is pointed
			%
			% A polyhedron { x | A*x<=b, C*x=d } is pointed if and only if
			% its lineality space null([A; C]) is empty.
			%
			% Literature:
			% http://www.ee.ucla.edu/ee236a/lectures/polyhedra.pdf
			
			answer = false(size(P));
			for i = 1:numel(P)
				% note that conversion to H-rep will be performed!
				answer(i) = isempty(null([P(i).A; P(i).Ae]));
			end
		end
		
		function R = and(P, Q)
			% P&Q computes intersection of two polytopes
			
			R = intersect(P, Q);
		end
		
		function [H, isConvex] = or(P, Q)
			% P|Q computes union of two polytopes
			%
			% U = P | Q returns a single polyhedron (the convex hull) if
			% the union of P and Q is convex. Otherwise U=[P; Q].
			%
			% [U, isConvex] = P | Q also returns a binary flag indicating
			% convexity.
			
			U = PolyUnion([P(:); Q(:)]);
			isConvex = U.isConvex();
			if isConvex
				H = U.convexHull();
			else
				H = [P(:); Q(:)];
			end
        end
        
        function [new, unique_idx, idx] = unique(P)
            % Returns unique components of an array of Polyhedron objects
            %
            % U = P.unique() returns only the unique elements of the array
            % P. [U, idx] = P.unique() returns the indices of unique
            % elements as the second output.
            
            idx = false(numel(P), numel(P));
            is_unique = true(size(P));
            for i = 1:numel(P)-1
                for j = i+1:numel(P)
                    if ~(idx(j, i))
                        answer = P(i)==P(j);
                        idx(j, i) = answer;
                        if answer
                            is_unique(j) = false;
                            break
                        end
                    end
                end
            end
            new = copy(P(is_unique));
            unique_idx = find(is_unique);
        end
        
        function Q = projectOnAffineHull(P, He)
            % Projects the polyhedron on its affine hull
            %
            % Given a polyhedron P = { x | A*x<=b, Ae*x=be }, calling
            % Q = P.projectOnAffineHull() produces a new polyhedron
            % Q = { z | A*F*z <= b - A*x0 } where F is in the null space of
            % Ae and x0 is any solution to Ae*x0=be.
            %
            % The equality constraints are identified via
            % Polyhedron/affineHull().
            %
            % The dimensionality of Q is dim(P)-rank(Ae, be).
            
            P.minHRep();
            Q = P.copy();
            for i = 1:numel(P)
                if nargin < 2
                    % to make Polyhedron/isBounded() more efficient, we
                    % allow to pass the pre-computed affine hull as the
                    % second argument 
                    He = P(i).affineHull();
                end
                if ~isempty(He)
                    % fully dimensional in a lower dimension
                    F = null(He(:, 1:end-1));
                    x0 = He(:, 1:end-1)\He(:, end);
                    A = P(i).A*F;
                    b = P(i).b - P(i).A*x0;
                    if isempty(A)
                        % projection is R^m with m = dim(P)-rank(Ae)
                        Q(i) = Polyhedron.fullSpace(size(A, 2));
                    else
                        Q(i) = Polyhedron(A, b);
                    end
                end
            end
                
        end
		
		function answer = isFullSpace(P)
			% returns true if the polyhedron represents R^n
			
			global MPTOPTIONS
			
			answer = false(size(P));
			for i = 1:numel(P)
				if P(i).hasHRep
					
					if ~isempty(P(i).He_int)
						% affine set is not R^n
					
					elseif isempty(P(i).He_int) && ...
						all(matNorm(P(i).H_int(:, 1:end-1)) < MPTOPTIONS.zero_tol) && ...
							all(P(i).H_int(:, end) > -MPTOPTIONS.zero_tol)
						% 0*x<=b with "b" non-negative => R^n
						%
						% Note that the polyhedron constructor already
						% removed rows like a'*x<=Inf, which are the only
						% other description of R^n
						answer(i) = true;
					end
					
				elseif P(i).hasVRep && ~isempty(P(i).R_int)
					% check whether rays span R^n
					
					% the rays span R^n if there exists lambda>=0 such
					% that R*lambda gives each basis vector of R^n
					
					answer(i) = true; % will be unset later if necessary
					
					% set of basis vectors
					E = [eye(P(i).Dim), -eye(P(i).Dim)];
					R = P(i).R_int';
					nR = size(R, 2);
					lp.f = zeros(nR, 1);
					lp.A = -eye(nR);
					lp.b = zeros(nR, 1);
					lp.Ae = R;
					lp.be = [];
					lp.lb = zeros(nR, 1);
					lp.ub = Inf(nR, 1);
					lp.quicklp = true;
					for j = 1:size(E, 2)
						% solve the feasibility LP:
						%   find lambda>=0, s.t. R*lambda=E(:, j)
						lp.be = E(:, j);
						res = mpt_solve(lp);
						if res.exitflag ~= MPTOPTIONS.OK
							% infeasible, the conic combinations of
							% rays failed to provide a basis vector
							answer(i) = false;
							break
						end
					end
					
				end
			end
        end
        
        function D = dual(obj)
            % Computes the polar dual of a polytope
            %
            % D = P.dual() computes the dual of the polyhedron P. The
            % polyhedron must contain the origin in its interior.
            %
            % The dual of P = { x | A1*x \le 1, A2*x \le 0 }
            % (which contains the origin in its interior) is 
            % D = { A1'y + A2*z | y \ge 0, \sum y=1, z \ge 0}, i.e.
            % the dual has the rows of A1 as its vertices (with the origin
            % included) and the rows of A2 as its rays.
            %
            % If the polyhedron is given in the vertex representation, i.e.,
            % P = { V'*y + R'*z | y \ge 0, \sum y = 1, z \ge 0 }, 
            % then the polar dual is D = { x | V*x \le 1, R*x \le 0 }.
            %
            % Special cases: the dual of an empty set is R^n and the dual
            % of R^n is the origin.
            %
            % Literature:
            % http://www.cis.upenn.edu/~cis610/polytope.ps
            
            global MPTOPTIONS
            if numel(obj)>1
                % element-wise operation on arrays
                D = obj.forEach(@(e) e.dual());
                return
            end
            
            if obj.isFullSpace()
                % polar dual of R^n is the origin
                D = Polyhedron(zeros(1, obj.Dim));
                
            elseif obj.isEmptySet()
                % polar dual of an empty set is the whole euclidian space
                D = Polyhedron.fullSpace(obj.Dim);
                
            else
                P = obj;
                origin = zeros(obj.Dim, 1);
                if ~obj.contains(origin)
                    error('The polyhedron must contain the origin in its interior.');
                end
                if obj.hasVRep
                    V = P.V;
                    D = Polyhedron([V; P.R], [ones(size(V, 1), 1); zeros(size(P.R, 1), 1)]);

                elseif obj.hasHRep
                    % split the H-rep to A1*x<=b1, A2*x<=0
                    A = obj.A; b = obj.b;
                    if ~isempty(obj.Ae)
                        % embed equalities as double-sided inequalities
                        A = [A; obj.Ae; -obj.Ae];
                        b = [b; obj.be; -obj.be];
                    end
                    iz = (abs(b)<=MPTOPTIONS.abs_tol);
                    A1 = A(~iz, :);
                    b1 = b(~iz);
                    A2 = A(iz, :);
                    % normalize A1*x<=b1 to A1*x<=1
                    if ~isempty(A1)
                        A1 = A1./repmat(b1, 1, obj.Dim);
                    end
                    A1 = [A1; origin']; % add the origin
                    
                    % compute the dual D = { A1'*y + A2'*z | y >= 0, sum(y)=1, z >= 0 }
                    D = Polyhedron('V', A1, 'R', A2);
                    
                else
                    % this should not happen, we should have at least one
                    % representation available
                    error('Unexpected error, please report to mpt@control.ee.ethz.ch');
                    
                end
            end
                
        end
        
        function x = randomPoint(obj)
            % Generates a random point inside a polytope
            %
            % Let v_i, i=1...n denote the vertices of a bounded polytope.
            % Then the random point is x = \sum L_i*v_i where 0<=L_i<=1 and
            % sum(L_i)=1.

            % TODO: more efficient procedure for H-polytopes
            % TODO: support polyhedra
            
            error(obj.rejectArray());
            if ~obj.isBounded()
                error('The polyhedron must be bounded.');
            end
            
            nv = size(obj.V, 1);
            % random vector with 0<=L(:)<=1 and sum(L)=1
            L = rand(nv, 1);
            L = L/sum(L);
            
            x = obj.V'*L;
        end
        
        function sol = fmax(obj, function_name)
            % Maximizes a function over a polyhedron
            %
            % sol = P.fmin(fun) maximizes a given function over a
            % polyhedron. The optimizer is returned in sol.xopt, the
            % objective value in sol.obj, and the optimization status in
            % sol.exitflag and sol.how.
            %
            % The function to be maximized must be attached to the
            % polyhedron and must be scalar-valued.
            %
            % If the problem to be solved is non-convex, YALMIP's status is
            % returned in sol.info.
            
            error(obj.rejectArray());
            
            if nargin < 2
                function_name = '';
            end
            sol = obj.fmin(function_name, -1);
        end
        
        function sol = fmin(obj, function_name, direction)
            % Minimizes a function over a polyhedron
            %
            % sol = P.fmin(fun) minimizes a given function over a
            % polyhedron. The optimizer is returned in sol.xopt, the
            % objective value in sol.obj, and the optimization status in
            % sol.exitflag and sol.how.
            %
            % The function to be minimized must be attached to the
            % polyhedron and must be scalar-valued.
            %
            % If the problem to be solved is non-convex, YALMIP's status is
            % returned in sol.info.
            
            global MPTOPTIONS

            error(obj.rejectArray());

            if nargin < 3
                direction = 1; % 1 = minimize, -1 = maximize
            end
            if nargin < 2
                function_name = '';
            end

            if isempty(function_name)
                % if no function is specified, take the first one
                fnames = obj.listFunctions();
                assert(numel(fnames)>0, 'The object has no functions.');
                assert(numel(fnames)==1, 'The object has multiple functions, specify the one to optimize.');
                function_name = fnames{1};
            end
            assert(ischar(function_name), 'The function name must be a string.');
            assert(obj.hasFunction(function_name), 'No such function "%s" in the object.', function_name);            

            fun = obj.Functions(function_name);

            if obj.isEmptySet()
                % trivially infeasible
                sol.xopt = NaN(obj.Dim, 1);
                sol.obj = direction*Inf;
                sol.exitflag = MPTOPTIONS.INFEASIBLE;
                sol.how = 'infeasible';
                return
            end
            
            if (isa(fun, 'QuadFunction') && direction==1 && min(eig(fun.H))>=0 ) || ...
                    (isa(fun, 'QuadFunction') && direction==-1 && max(eig(fun.H))<=0 ) || ...
                    isa(fun, 'AffFunction') || ...
                    ((isa(fun, 'OneNormFunction') || isa(fun, 'InfNormFunction')) && direction==1)
                % min of convex quadratic function => QP
                % max of concave quadratic function => Qp
                % min/max of linear function => LP
                % min of 1-norm or Inf-norm => LP
                
                assert(fun.R==1, 'The function to minimize must be scalar.');
                
                prob.lb = [];
                prob.ub = [];
                if isa(fun, 'QuadFunction')
                    prob.H = direction*2*fun.H; % mpt_solve minimizes 0.5 x'Hx + f'x !
                    prob.quickqp = true;
                else
                    prob.quicklp = true;
                end
                if isa(fun, 'OneNormFunction')
                    % min ||Qx||_1 is equivalent to
                    % min sum(e) s.t. -e <= Qx <= e
                    % where "e" is a rows(Q)-by-1 vector
                    Q = fun.weight;
                    ne = size(Q, 1);
                    prob.A = [obj.A, zeros(size(obj.A, 1), ne); ...
                        Q, -eye(ne); ...
                        -Q, -eye(ne)];
                    prob.b = [obj.b; zeros(2*ne, 1)];
                    prob.Ae = [obj.Ae, zeros(size(obj.Ae, 1), ne)];
                    prob.be = obj.be;
                    prob.f = [zeros(1, obj.Dim), ones(1, ne)];
                    
                elseif isa(fun, 'InfNormFunction')
                    % min ||Qx||_inf is equivalent to
                    % min e s.t. -e <= Qx <= e
                    % where "e" is a scalar
                    Q = fun.weight;
                    nQ = size(Q, 1);
                    prob.A = [obj.A, zeros(size(obj.A, 1), 1); ...
                        Q, -ones(nQ, 1); ...
                        -Q, -ones(nQ, 1)];
                    prob.b = [obj.b; zeros(2*nQ, 1)];
                    prob.Ae = [obj.Ae, zeros(size(obj.Ae, 1), 1)];
                    prob.be = obj.be;
                    prob.f = [zeros(1, obj.Dim), 1];
                    
                else
                    % add linear terms for LPs and QPs
                    prob.f = direction*fun.F;
                    prob.A = obj.A;
                    prob.b = obj.b;
                    prob.Ae = obj.Ae;
                    prob.be = obj.be;
                    
                end
                
                sol = mpt_solve(prob);
                switch sol.exitflag
                    case MPTOPTIONS.OK
                        if isa(fun, 'NormFunction')
                            % remove the slacks from the optimizer
                            sol.xopt = sol.xopt(1:obj.Dim);
                        else
                            sol.obj = sol.obj + direction*fun.g;
                        end
                        
                    case MPTOPTIONS.INFEASIBLE
                        sol.obj = Inf;
                        sol.xopt = NaN(obj.Dim, 1);
                        
                    case MPTOPTIONS.UNBOUNDED
                        sol.obj = -Inf;
                        sol.xopt = NaN(obj.Dim, 1);
                        
                    otherwise
                        error('Unknown result from the solver.');
                end
                
            else
                % general functions are minimized via YALMIP
                
                x = sdpvar(obj.Dim, 1);
                con = [obj.A*x<=obj.b, obj.Ae*x<=obj.be];
                J = direction*fun.feval(x);
                assert(isscalar(J), 'The function to minimize must be scalar.');
                info = solvesdp(con, J, sdpsettings('verbose', 0, 'debug', 1));
                switch info.problem
                    case {0, 3, 4, 5}
                        sol.xopt = double(x);
                        sol.obj = double(J);
                        sol.how = 'ok';
                        sol.exitflag = MPTOPTIONS.OK;
                        
                    case {2, 12, 15}
                        % unbounded (note that we know the problem should
                        % have been feasible since the domain is not an
                        % empty set)
                        sol.obj = -Inf;
                        sol.xopt = NaN(obj.Dim, 1);
                        sol.how = 'unbounded';
                        sol.exitflag = MPTOPTIONS.UNBOUNDED;
                        
                    otherwise
                        % infeasible
                        sol.obj = Inf;
                        sol.xopt = NaN(obj.Dim, 1);
                        sol.how = 'infeasible';
                        sol.exitflag = MPTOPTIONS.INFEASIBLE;
                end
                
                sol.info = info;
            end
            
            % adjust the cost when maximizing
            sol.obj = sol.obj*direction;
        end
        
        function [b, bi, wi] = barycentricCoordinates(P)
            % Computes barycentric coordinate functions for a polytope P
            %
            % Given a polytope P = convh({v1, ..., vm}), this function computes the
            % barycentric coordinate functions bi(x), i=1,...,m such that:
            %   bi(x) >= 0           (non-negativity)
            %   \sum_i bi(x) = 1     (partition of unity)
            %   \sum vi*bi(x) = vi   (linear precision)
            %
            % Usage: [b, bi, wi] = barycentricCoordinates(P)
            %
            % Outputs:
            %    b: function handle that takes a point "x" as an input, and
            %       returns a vector of values bi(x)
            %   bi: cell array of barycentric functions, one for each
            %       vertex (each function takes "x" as an input)
            %   wi: cell array of weighting functions, one for each vertex
            %
            % Limitation: each vertex of P must be simple, i.e., it must have exactly D
            % incident facets where D is the dimension of P.
            %
            % Literature:
            %   J. Warren, S. Schaefer, A. Hirani, and M. Desbrun:
            %   Barycentric coordinates for convex sets.
            %   Advances in Computational Mathematics, 27(3):319-338, Aug 2007.
            
            global MPTOPTIONS
            
            assert(numel(P)==1, 'Single polytope please.');
            assert(P.isBounded(), 'The polytope must be bounded');
            
            P.normalize();
            imap = P.incidenceMap();
            V = imap.V';
            nv = size(V, 2);
            nx = P.Dim;
            % shift to origin, enlarge and shift back to guarantee
            % partition of unity at the vertices
            xc = repmat(P.chebyCenter().x, 1, nv);
            V = (1+MPTOPTIONS.rel_tol)*(V - xc) + xc;
            
            % Compute weight functions wi(x), see Section 2 of the paper
            wi = {};
            for i = 1:nv
                % indices of incident facets
                v = V(:, i);
                ind_v = find(imap.incVH(i, :));
                assert(length(ind_v)==nx, sprintf('Vertex no. %d is not simple.', i));
                % facets incident to the vertex
                nind_v = P.A(ind_v, :);
                % volume of the parallelepiped span by the normals to the facets
                % incident to the vertex
                kappa = abs(det(nind_v));
                % numerator of (4)
                z = '@(x) kappa/((nind_v(1, :)*(v-x))';
                for j = 2:size(nind_v, 2)
                    z = [z, sprintf('*(nind_v(%d, :)*(v-x))', j)];
                end
                z = [z, ')'];
                wi{i} = eval(z);
            end
            
            bi = cell(1, nv);
            for i = 1:nv
                bi{i} = @(x) (wi{i}(x))/sum(cellfun(@(e) e(x), wi));
            end
            b = @(x) cellfun(@(e) e(x), bi);
        end
        
        function gi = interpolate(P, g)
            % Interpolates values or functions at vertices over a polytope
            %
            % If G is a matrix with "nv" columns where "nv" is the number
            % of vertices of a polytope P, then Gi=interpolate(P, G)
            % returns an interpolation function Gi(X) that interpolates
            % the values at vertices (represented by the columns of G) at
            % a point X of the polytope P using its barycentric
            % coordinates.
            %
            % If G is a cell array of function handles, Gi(X) is a function
            % which interpolates the values of G at individual vertices at
            % a point X of P.
            %
            % Example:
            %   P = Polyhedron.unitBox(2); % box with 4 vertices
            %   G = [1, 2, 3, 4];          % values at vertices
            %   Gi = P.interpolate(G);     % obtain the interpolation
            %   x = [0.2; 0.6];            % point in P
            %   Gix = Gi(x)                % 2.36 as the interpolated value
            %
            % Interpolation with functions defined over vertices:
            %   G = { @(z) z+1, @(z) sin(z), @(z) z^2, @(z) z };
            %   Gi = P.interpolate(G);
            %   z0 = 0.4;
            %   x = [0.2; 0.6];
            %   Gix = Gi(x, z0)   % 0.4381 as the interpolated value
            %
            % See also Polyhedron/barycentricCoordinates
            
            assert(numel(P)==1, 'Single polyhedron please.');
            
            function gi = baryinterp(x, barycoords, g, varargin)
                f = cellfun(@(e) e(varargin{:}), g, 'UniformOutput', false);
                gi = sum(barycoords(x).*cat(2, f{:}), 2);
            end

            barycoords = P.barycentricCoordinates();
            nv = size(P.V, 1);
            if iscell(g)
                % interpolate functions
                assert(numel(g)==nv, 'The cell array must have %d elements.', nv);
                na = nargin(g{1});
                if na==0
                    % functions with no inputs
                    gi = @(x) baryinterp(x, barycoords, g);
                else
                    % functions with inputs
                    args = 'z1';
                    for i = 2:na
                        args = [args, sprintf(', z%d', i)];
                    end
                    gi = eval(sprintf('@(x, %s) baryinterp(x, barycoords, g, %s)', args, args));
                end
            elseif isnumeric(g)
                % interpolate values at vertices
                assert(size(g, 2)==nv, 'The matrix must have %d columns.', nv);
                gi = @(x) sum(barycoords(x).*g, 2);
            else
                error('The second input can be only a cell array of function handles, or a matrix.');
            end
        end
        
	end
	
	methods(Hidden)
		% non-public methods

		% Polyhedron/isInside() implicitly assumes that the
		% H-representation is available. Use it only if you know what you
		% are doing! In general cases you should use Polyhedron/contains()
		% instead.
		[isin, inwhich, closest] = isInside(P, x, Options);

        function answer = doesIntersect(P1, P2, how)
            % Returns true if P1 and P2 intersect
            %
            % answer = P1.doesIntersect(P2) returns true if the intersection
            % of P1 and P2 is not an empty set.
            %
            % answer = P1.doesIntersect(P2, 'fully') returns true if the
            % intersection is a fully dimensional set.
            %
            % P1 and P2 must be single Polyhedron objects. The check
            % requires the H-representation of both polyhedrons.
            %
            % This is a faster version of P1.intersect(P2).isEmptySet() and
            % P1.intersect(P2).isFulDim().
            
            % Note that this method is intended to be as fast as possible
            % for a good performance of PolyUnion/min()
            
            global MPTOPTIONS
            if nargin<3
                how = '';
            end
            
            assert(isa(P1, 'Polyhedron') && isa(P2, 'Polyhedron'), 'Both objects must be a Polyhedron.');
            assert(numel(P1)==1 && numel(P2)==1, 'Arrays are not supported.');
            assert(P1.Dim==P2.Dim, 'Both polyhedra must be in the same dimension.');
            
            % first check intersection of bounding boxes
            if ~(isfield(P1.Internal, 'lb') && isfield(P1.Internal, 'ub'))
                P1.outerApprox();
            end
            if ~(isfield(P2.Internal, 'lb') && isfield(P2.Internal, 'ub'))
                P2.outerApprox();
            end
            lb1 = P1.Internal.lb;
            ub1 = P1.Internal.ub;
            lb2 = P2.Internal.lb;
            ub2 = P2.Internal.ub;
            if any(lb1-ub2>10*MPTOPTIONS.rel_tol) || ...
                    any(lb2-ub1>10*MPTOPTIONS.rel_tol)
                % bonding boxes do not intersect => sets do not intersect
                answer = false;
                return
            end
            
            P1.computeHRep();
            P2.computeHRep();
            
            % check emptiness/full dimensionality of the intersection
            if isempty([P1.H_int; P1.He_int]) || isempty([P2.H_int; P2.He_int])
                % one of them is an empty set => no intersection
                answer = false;
                return
            end
            
            H = [P1.H_int; P2.H_int; P1.He_int; -P1.He_int; P2.He_int; -P2.He_int];
            if isequal(how, 'fully')
                answer = fast_isFullDim(H);
            else
                answer = ~fast_isEmptySet(H);
            end
        end
        
	end
	
	methods(Static)
	
		function S = emptySet(dim)
			% Polyhedron.emptySet(n) constructs an empty set in R^n

			narginchk(1, 1);
			S = Polyhedron(zeros(0, dim), zeros(0, 1));
			S.Internal.Empty = true;
			S.Internal.FullDim = false;
			S.Internal.lb = Inf(dim, 1);
			S.Internal.ub = -Inf(dim, 1);
		end

		function S = fullSpace(dim)
			% Polyhedron.fullSpace(n) constructs the H-representation of R^n

			% R^n is represented by 0'*x<=1
			S = Polyhedron(zeros(1, dim), 1);
			S.irredundantHRep = true;
			S.Internal.Empty = (dim==0); % R^0 is an empty set
			S.Internal.FullDim = (dim>0); % R^0 is not fully dimensional
			S.Internal.Bounded = (dim==0); % R^0 is bounded
			S.Internal.lb = -Inf(dim, 1);
			S.Internal.ub = Inf(dim, 1);
		end
		
		function B = unitBox(dim)
			% Polyhedron.unitBox(n) constructs a unit box in "n" dimensions
			
			narginchk(1, 1);
			B = Polyhedron('lb', -ones(dim, 1), 'ub', ones(dim, 1), ...
				'irredundantHRep', true);
		end

		function new = loadobj(obj)
			% load Polyhedron objects
			
			% NOTE: loadobj() must be a static method!
			
			% empty polyhedron will provide us with the version in
			% base.version
			base = Polyhedron;
			if ~isprop(obj, 'version')
				% ancient version before September 25, 2012
				if isempty(obj.H) && isempty(obj.He)
					% Vrep
					new = Polyhedron('V', obj.V, 'R', obj.R);
				elseif isempty(obj.V) && isempty(obj.R)
					% Hrep
					new = Polyhedron('H', obj.H, 'He', obj.He);
				else
					% mixed V/H rep
					new = Polyhedron('H', obj.H, 'He', obj.He, ...
						'V', obj.V, 'R', obj.R);
				end
				% TODO: import functions?
			elseif obj.version > base.version
				error('Version of loaded object: %.1f, supported version: %.1f.', ...
					obj.version, base.version);
			else
				% object is up-to-date, just use the copy constructor
				new = Polyhedron(obj);
			end
		end
	end
	
	methods (Access=protected)
		
		% function prototypes (plotting methods must remain protected)
		h = fplot_internal(obj, function_name, options)
		h = plot_internal(obj, options)

	end
	
end

