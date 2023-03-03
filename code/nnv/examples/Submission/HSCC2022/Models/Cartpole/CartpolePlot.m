function [] = CartpolePlot(allReach, init_state)
%CARTPOLEPLOT Plot the results of the reachability analysis

    f1 = figure;
    hold on;
    grid;
    set(gcf,'Color',[1 1 1]);
    set(gca, 'GridAlpha', 1); % Set transparency of grid
    %set(gca, 'color', [17 17 17]/19); % Set background color 
    set(gcf,'inverthardcopy','off'); % Enable saving the figure as it is
    Star.plotBoxes_2D_noFill(allReach.step_sets,3,4,'b');
    Star.plotBoxes_2D_noFill(allReach.int_reachSet,3,4,'b');
    xlabel('$\theta$ (in rad)','Interpreter','latex');
    ylabel('$\dot{\theta}$ (in rad/s)','Interpreter','latex');
    title('Phase portrait (angle)','Interpreter','latex');
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

    f2 = figure;
    hold on;
    grid;
    set(gcf,'Color',[1 1 1]);
    set(gca, 'GridAlpha', 1); % Set transparency of grid
    %set(gca, 'color', [17 17 17]/19); % Set background color 
    set(gcf,'inverthardcopy','off'); % Enable saving the figure as it is
    Star.plotBoxes_2D_noFill(allReach.step_sets,1,2,'b');
    Star.plotBoxes_2D_noFill(allReach.int_reachSet,1,2,'b');
    xlabel('$x$ (in m)','Interpreter','latex');
    ylabel('$\dot{x}$ (in m/s)','Interpreter','latex');
    title('Phase portrait (position)','Interpreter','latex');
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

