classdef Conversion
    %Conversion class contains some basic conversion methodd for merging polyhedra 
    %   Dung Tran
    
    properties
    end
    
    methods(Static)
         % convert equalities to inequalities for a polyhedron
         function new_P = removeEqualities(P)
             % @P: a Polyhedron
             % @new_P: a new polyhedron without equalities
                         
             % author: Dung Tran
             % date: 11/5/2018
             
             if ~isa(P, 'Polyhedron')
                 error('Input set is not a polyhedron');
             end
             
             if isempty(P.Ae)
                 new_P = P;
             else
                 A = vertcat(P.A, P.Ae, -P.Ae);
                 b = vertcat(P.b, P.be, -P.be);
                 new_P = Polyhedron('A', A, 'b', b);
             end

         end
         
         % concatenate two polyhedra to make a higher dimensional
         % polyhedron
         function P = concatenatePolyhedron(Ps)
             % @Ps: an array of input polyhedron
             % @P: the concatenated polyhdron P: [x; y], x \in P1, y \in P2
             
             % author: Dung Tran
             % date: 11/5/2018
             
             
             n = length(Ps);
             A = [];
             b = [];
             for i=1:n
                 P1 = Conversion.removeEqualities(Ps(i));
                 A = blkdiag(A, P1.A);
                 b = vertcat(b, P1.b);                
             end
             
             P = Polyhedron('A', A, 'b', b);
             
         end
         
         % to Star
         function S = toStar(P)
             % convert a Polyhedron to a Star
             % Author: Dung Tran
             % Date: 11/16/2018
             
             if ~isa(P, 'Polyhedron')
                 error('Input is not a polyhedron');
             end
             
             dim = P.Dim;
             P1 = Conversion.removeEqualities(P);
             
             c = zeros(dim, 1);
             V = eye(dim);
             
             S = Star([c V], P1.A, P1.b);
             nV = S.nVar; 
             S.predicate_lb = -ones(nV,1);
             S.predicate_ub = ones(nV,1);             
         end
         
         % to Box
         function B = toBox(P)
             % convert a Polyhedron to a Box
             if ~isa(P, 'Polyhedron')
                 error('Input is not a polyhedron');
             end
             
             P.outerApprox;
             lb = P.Internal.lb;
             ub = P.Internal.ub;
             B = Box(lb, ub);

         end
         
         % get range of Polyhedron at specific index
         function [xmin, xmax] = getRange(Ps, index)
             % @P: an array of polyhedra.
             % @index: index
             
             % Dung Tran
             % date: 11/16/2018
             
             n = length(Ps); % number of polyhedra
             xmin = zeros(n,1);
             xmax = zeros(n,1);
             for i=1:n 
                 
                P = Ps(i);
                if ~isa(P, 'Polyhedron')
                   error('Input is not a polyhedron');
                end
    
                dim = P.Dim;
                f = zeros(1, dim);
                f(index) = 1;
                options = optimset('Display','none');
                [~, fval, exitflag, ~] = linprog(f, P.A, P.b, P.Ae, P.be, [], [], [], options);
            
                if exitflag > 0
                    xmin(i) = fval;
                else
                    error('Cannot find an optimal solution');
                end          

                [~, fval, exitflag, ~] = linprog(-f, P.A, P.b, P.Ae, P.be, [], [], [], options);
                if exitflag > 0
                    xmax(i) = -fval;
                else
                    error('Cannot find an optimal solution');
                end
                
             end
             
                          
         end
         
         
    end
    
end

