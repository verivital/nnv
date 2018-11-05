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
         function P = concatenatePolyhedron(P1, P2)
             % @P1: the first polyhedron
             % @P2: the second polyhedron
             % @P: the concatenated polyhdron P: [x; y], x \in P1, y \in P2
             
             % author: Dung Tran
             % date: 11/5/2018
             
             newP1 = Conversion.removeEqualities(P1);
             newP2 = Conversion.removeEqualities(P2);
             
             [nP1, mP1] = size(newP1.A);
             [nP2, mP2] = size(newP2.A);
             
             Z1 = zeros(nP1, mP2);
             Z2 = zeros(nP2, mP1);
             
             A = [newP1.A Z1; Z2 newP2.A];
             b = [newP1.b; newP2.b];
             
             P = Polyhedron('A', A, 'b', b);
                          
         end
    end
    
end

