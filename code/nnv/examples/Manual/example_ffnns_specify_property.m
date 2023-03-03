
%/* Example of specifying a property of FFNN */

% unsafe region: y1 >= 5 (the network has two outputs)

G = [-1 0]; % condition matrix

g = -5;     % condition vector

U = HalfSpace(G, g); % unsafe region object