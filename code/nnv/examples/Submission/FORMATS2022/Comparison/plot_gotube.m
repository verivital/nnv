function plot_gotube(sets,x1,x2,dim,clr)
%PLOT_GOTUBE reach sets
hold on;
% Use viscircles to plot Gotube reach sets
centers = [sets(x1+1,:);sets(x2+1,:)]';
radii = sets(dim+2,:)';
viscircles(centers, radii,'Color',clr);
end

