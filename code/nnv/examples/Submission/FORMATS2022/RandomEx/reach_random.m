function reach_random(neuralode,x0,fname,unc)
dim = length(x0);
% Create initial state
lb = x0-unc;
ub = x0+unc;
R0 = Star(lb,ub);

y = neuralode.evaluate(x0); % Simulation
t = tic;
R = neuralode.reach(R0); % Reachability
rT = toc(t);
rname = "results/"+string(fname)+"_"+string(unc)+".mat";
save(rname,'R','rT');


if dim > 1
    plot_combos = nchoosek(1:dim,2);
else
    plot_combos = [1 1];
end

% Plot results
for i=1:size(plot_combos,1)
    tmp = plot_combos(i,:);
    ffname = "figures/"+string(fname)+"_"+string(unc)+"_"+string(tmp(1))+string(tmp(2))+".pdf";
    f = figure;
    hold on;
    Star.plotBoxes_2D_noFill(R,tmp(1),tmp(1),'k');
    pg = plot(y(tmp(1),1),y(tmp(1),1),'k');
    xlabel("x_"+string(tmp(1)));
    ylabel("x_"+string(tmp(2)));
    ax = gca; % Get current axis
    ax.XAxis.FontSize = 15; % Set font size of axis
    ax.YAxis.FontSize = 15;
    legend(pg,{'NNVODE (ours)'},"Location","best",'FontSize',14);
    exportgraphics(f,ffname,'ContentType','vector');
end
end

