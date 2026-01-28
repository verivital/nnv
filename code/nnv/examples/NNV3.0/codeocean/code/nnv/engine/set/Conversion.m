classdef Conversion
    %Conversion class contains some basic conversion methodd for merging polyhedra 
    %   Dung Tran
    
    properties % none

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
         
         % concatenate two polyhedra to make a higher dimensional polyhedron
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
         
         % Polyhedron to Star
         function S = toStar(P)
             % convert a Polyhedron to a Star
             % Author: Dung Tran
             % Date: 11/16/2018
             
             if ~isa(P, 'Polyhedron')
                 error('Input is not a polyhedron');
             end
             
             dim = P.Dim;
             P1 = Conversion.removeEqualities(P);
             c = cast(zeros(dim, 1), 'like', P1.A);
             V = cast(eye(dim), 'like', P1.A);
             P1.outerApprox;
             pred_lb = P1.Internal.lb;
             pred_ub = P1.Internal.ub;  
             S = Star([c V], P1.A, P1.b, pred_lb, pred_ub);
         end
         
         % Polyhedron to Box
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

             if ~isa(Ps, 'Polyhedron')
               error('Input is not a polyhedron');
            end
             
             n = length(Ps); % number of polyhedra
             xmin = cast(zeros(n,1), 'like', Ps(1).A);
             xmax = cast(zeros(n,1), 'like', Ps(1).A);
             for i=1:n 
                P = Ps(i);
                dim = P.Dim;
                f = cast(zeros(1, dim), 'like', P.A);
                f(index) = 1;
                
                [fval, exitflag] = lpsolver(f,P.A, P.b, P.Ae, P.be, [], []);
                if strcmp(exitflag, "l1")
                    xmin(i) = fval;
                else
                    error('Cannot find an optimal solution');
                end          

                [fval, exitflag] = lpsolver(-f, P.A, P.b, P.Ae, P.be, [], []);
                if strcmp(exitflag, "l1")
                    xmax(i) = -fval;
                else
                    error('Cannot find an optimal solution');
                end
                
             end
                          
         end
         
    end
    
end

