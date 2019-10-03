load reachSet_polyhedron.mat;
R1 = S; % polyhedron reach sets (computed by polyhedron method)
load reachSet_Interval.mat;
R2 = S; % Interval reach sets (computed by interval method)
load reachSet_Star.mat;
R3 = X_cell; % exact star reach sets (computed by exact star set method)
load reachSet_Star_approx.mat;
R4 = S; % approximate star reach sets (computed by approximate star set method)


n1 = length(R1);
n2 = length(R2);
n3 = length(R3);
n4 = length(R4); 

map = [1 0 0; 0 1 0];

R1_map = cell(1, n1);
for i=1:n1
    Ri = R1{1,i};
    m = length(Ri);
    Y = [];
    for j=1:m
        Y = [Y Ri(j).affineMap(map)];
    end
    R1_map{1, i} = Y;
end

R2_map = cell(1, n2);
for i=1:n2
    Ri = R2{1,i};
    m = length(Ri);
    Y = [];
    for j=1:m
        Y = [Y Ri(j).affineMap(map)];
    end
    R2_map{1, i} = Y;
end


R3_map = cell(1, n3);
for i=1:n3
    Ri = R3{1,i};
    m = length(Ri);
    Y = [];
    for j=1:m
        Y = [Y Ri(j).affineMap(map, [])];
    end
    R3_map{1, i} = Y;
end

R4_map = cell(1, n4);
for i=1:n4
    Ri = R4(i);
    m = length(Ri);
    Y = [];
    for j=1:m
        Y = [Y Ri(j).affineMap(map, [])];
    end
    R4_map{1, i} = Y;
end





% compute bounds of polyhedron reach sets
lb_poly = [];
ub_poly = [];
for i=1:n1
    B = Reduction.hypercubeHull(R1_map{1,i});
    lb_poly = [lb_poly B.lb];
    ub_poly = [ub_poly B.ub];
end

% compute bounds of interval reachable sets
lb_interval = [];
ub_interval = [];
for i=1:n2
    B = Reduction.hypercubeHull(R2_map{1, i});
    lb_interval = [lb_interval B.lb];
    ub_interval = [ub_interval B.ub];
end

% compute exact bounds
lb_star_exact = [];
ub_star_exact = [];
for i=1:n3
    B = Star.get_hypercube_hull(R3_map{1,i});
    lb_star_exact = [lb_star_exact B.lb];
    ub_star_exact = [ub_star_exact B.ub];
end

% compute exact bounds
lb_star_approx = [];
ub_star_approx = [];
for i=1:n4
    B = Star.get_hypercube_hull(R4_map{1,i});
    lb_star_approx = [lb_star_approx B.lb];
    ub_star_approx = [ub_star_approx B.ub];
end

save errors.mat lb_poly ub_poly lb_interval ub_interval lb_star_exact ub_star_exact lb_star_approx ub_star_approx 





%subplot(2,2,1);
%for i=1:n1
%    Si = R1_map{1, i};
%    plot(Si, 'Color', 'red');
%    hold on;
%end
%xlabel('Distance (m)');
%ylabel('Velocity (m/s)');
%title('Polyhedron');
%set(gca,'FontSize',14)

