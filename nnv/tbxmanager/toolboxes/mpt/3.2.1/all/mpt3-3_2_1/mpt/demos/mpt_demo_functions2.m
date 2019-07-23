function mpt_demo_functions2
%
% mpt_functions_demo2
%

%% demo that illustrates functions over unions of  polyhedra
close all


%% constructing union of triangular polyhedra
disp('Create random polyhedron')
P = 10*ExamplePoly.randVrep

disp(' '); disp(' ');
disp('Triangulate the polyhedron to get a complex.');
T = P.triangulate

disp(' '); disp(' ');
disp('For each of the polyhedron, assign affine function with the same name')
for i=1:numel(T)
    T(i).addFunction(AffFunction(eye(2),[-1;1]),'phi');
end

pause
disp(' '); disp(' ');
disp('Construct the polyunion object U.')
U = PolyUnion('Set',T,'FullDim',true,'Bounded',true,'Overlaps',false,'Convex',true)

pause
disp(' '); disp(' ');
disp('Plot the function over the polyhedra')
U.fplot

pause,
close all

%% construct overlapping union
disp('Create 3 random polyhedra.');
for i=1:3
    Q(i) = ExamplePoly.randVrep+5*rand(2,1);
end

pause
disp('Assign two quadratic functions to each of the polyhedra.')
for i=1:3
    Q(i).addFunction(QuadFunction(eye(2),randn(1,2),randn(1)),'alpha')
    Q(i).addFunction(QuadFunction(eye(2),randn(1,2),randn(1)),'beta')
end
pause

disp(' '); disp(' ');
disp('Create union without specifying any properties.');
PU = PolyUnion(Q)

pause
disp(' '); disp(' ');
disp('Plot the functions over polyhedra.');
PU.fplot('beta')

pause
disp('Plot the functions over polyhedra with some properties.');
PU.fplot('beta','show_set',true)


end
