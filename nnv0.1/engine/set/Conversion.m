classdef Conversion
    %Conversion class contains some basic conversion method for merging polyhedra 
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
    end
    
end

