function [] = ACASXuPlot(allReach, init_state)
%ACASXUPLOT Plot the results of the reachability analysis

    function plot_relative_ACASXu(stars, color)
    
            n = length(stars);
            x_own_min = zeros(n,1);
            x_own_max = zeros(n,1);
            y_own_min = zeros(n,1);
            y_own_max = zeros(n,1);
            x_int_min = zeros(n,1);
            x_int_max = zeros(n,1);
            y_int_min = zeros(n,1);
            y_int_max = zeros(n,1);
    
            for i=1:n
                if isa(stars(i), 'Star')
                    [x_own_min(i), x_own_max(i)] = stars(i).getRange(1);
                    [y_own_min(i), y_own_max(i)] = stars(i).getRange(2);
                    [x_int_min(i), x_int_max(i)] = stars(i).getRange(4);
                    [y_int_min(i), y_int_max(i)] = stars(i).getRange(5);
                else
                    error('%d th input object is not a star', i);
                end
            end
            
            xmin = zeros(n,1);
            xmax = zeros(n,1);
            ymin = zeros(n,1);
            ymax = zeros(n,1);
            
            for i=1:n
    
                xmin(i) = x_int_min(i) - x_own_max(i);
                xmax(i) = x_int_max(i) - x_own_min(i);
                ymin(i) = y_int_min(i) - y_own_max(i);
                ymax(i) = y_int_max(i) - y_own_min(i);
    
            end
    
    
            for i=1:n
    
                x = [xmin(i) xmax(i) xmax(i) xmin(i) xmin(i)];
                y = [ymin(i) ymin(i) ymax(i) ymax(i) ymin(i)];
    
                hold on;
                plot(x, y, color);
    
            end
    
    
    end

    f = figure;
    hold on;
    grid;
    set(gcf,'Color',[1 1 1]);
    set(gca, 'GridAlpha', 1); % Set transparency of grid 
    set(gcf,'inverthardcopy','off'); % Enable saving the figure as it is
    plot_relative_ACASXu(allReach.step_sets,'b');
    plot_relative_ACASXu(allReach.int_reachSet,'b');
    ownship = plot(0,0,'g*','linewidth', 4);
    intruder = plot(init_state(4),init_state(5),'b*','linewidth', 4);
    legend([ownship,intruder],{'ownship','intruder'});
    xlabel('$x$ (in ft)','Interpreter','latex');
    ylabel('$y$ (in ft)','Interpreter','latex');
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

