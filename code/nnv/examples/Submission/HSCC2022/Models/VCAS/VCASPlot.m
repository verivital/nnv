function [] = VCASPlot(allReach, init_state)
%VCASPLOT Plot the results of the reachability analysis

    f = figure;
    hold on;
    grid;
    set(gcf,'Color',[1 1 1]);
    set(gca, 'GridAlpha', 1); % Set transparency of grid
    set(gcf,'inverthardcopy','off'); % Enable saving the figure as it is
    Star.plotBoxes_2D_noFill(allReach.step_sets,4,1,'b');
    Star.plotBoxes_2D_noFill(allReach.int_reachSet,4,1,'b');
    ownship = plot(0,0,'g*','linewidth', 4);
    intruder = plot(init_state(4), init_state(1),'b*','linewidth', 4);
    legend([ownship,intruder],{'ownship','intruder'});
    xlabel('$\tau$ time before loss of horizontal separation in s','Interpreter','latex');
    ylabel('$h$ relative altitude in ft','Interpreter','latex');
    title('Trajectory of intruder relative to ownship','Interpreter','latex');
    ax = gca; % Get current axis
    ax.GridColor = 'w'; % Set grid lines color
    ax.XAxis.FontSize = 15; % Set font size of axis
    ax.YAxis.FontSize = 15;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
end

