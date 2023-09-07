%% Show rotation examples

% Let's define a 2-dimensional shape (rectangle)

% center at [0.5, 0.5]
% side lengths = 1
V = [0.5 0.5 0; 0.5 0 0.5];
C = [0,0];
d = 0;
pred_lb = [-1;-1];
pred_ub = [1;1];

X = Star(V,C,d,pred_lb, pred_ub);

figure;
Star.plots(X,'b');

% Now, if we want to rotate this set 45 degrees, what would it look like?

% the center needs to be the same
V = [0, sqrt(2)/4, -sqrt(2)/4;
sqrt(2)/2, sqrt(2)/4, sqrt(2)/4];
C = [0 0];
d = [0];
pred_lb = [-1;-1];
pred_ub = [ 1; 1];

X2 = Star(V,C,d,pred_lb, pred_ub);

figure;
Star.plots(X2,'b');
% xlim([-3 3])
% ylim([-3 3])


%% Let's try the rotation here as well

% rotate the rectangle by an angle
angle = 45;
angle = deg2rad(angle);

rot_mat = [cos(angle), -sin(angle);
    sin(angle), cos(angle)];

X_r = X.affineMap(rot_mat, []);

figure;
Star.plots(X_r,'b');

% X_r = X2

% We can see that for a simple set, this actually works
% so, why doesn't it work in general?

% This is the function to convert a Star to polyhedron
% We can see that the lower and upper bounds are just added constraints, so
% we should be able to simply represent the Star using C and d and remove
% the predicate bounds
%
% function P = toPolyhedron(obj)
%     
%     b = obj.V(:, 1);        
%     W = obj.V(:, 2:obj.nVar + 1);
%     
%     if ~isempty(obj.predicate_ub)
%         C1 = cast([eye(obj.nVar); -eye(obj.nVar)], 'like', obj.V);
%         d1 = [obj.predicate_ub; -obj.predicate_lb];
%         Pa = Polyhedron('A', [obj.C;C1], 'b', [obj.d;d1]);
%         P = W*Pa + b;
%     else
%         Pa = Polyhedron('A', [obj.C], 'b', [obj.d]);
%         P = W*Pa + b;
%     end
% end

%% Let's add constraints to the C and d variables rather than the predicate lower and uper bounds

% the center needs to be the same
V = [0, sqrt(2)/4, -sqrt(2)/4;
sqrt(2)/2, sqrt(2)/4, sqrt(2)/4];
C = [0 0];
d = [0];
pred_lb = [-1;-1];
pred_ub = [ 1; 1];

X2 = Star(V,C,d,pred_lb, pred_ub);