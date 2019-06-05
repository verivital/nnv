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


figure; 

subplot(2,2,1);
for i=1:n1
    Si = R1_map{1, i};
    plot(Si, 'Color', 'red');
    hold on;
end
xlabel('Distance (m)');
ylabel('Velocity (m/s)');
title('Polyhedron');
set(gca,'FontSize',14)

subplot(2,2,2);
for i=1:n2
    Si = R2_map{1, i};
    Si.plot;
    hold on;
end
xlabel('Distance (m)');
ylabel('Velocity (m/s)');
title('Interval');
set(gca,'FontSize',14)

subplot(2,2,3);
for i=1:n3
    Si = R3_map{1, i};
    Star.plots(Si);
    hold on;
end
xlabel('Distance (m)');
ylabel('Velocity (m/s)');
title('Exact-Star');
set(gca,'FontSize',14)


subplot(2,2,4);
for i=1:n4
    Si = R4_map{1, i};
    Star.plots(Si);
    hold on;
end
xlabel('Distance (m)');
ylabel('Velocity (m/s)');
title('Approximate-Star');
set(gca,'FontSize',14)
