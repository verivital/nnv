function create_figure()

    % Load data
    cd N_2_0;
    n2 = load("N2_0_results.mat");
    cd ..;
    
    cd N_4_4;
    n4 = load("N4_4_results.mat");
    cd ..;
    
    cd N_8_0;
    n8 = load("N8_0_results.mat");
    cd ..;
    
    % Process results
    time_n2 = sum(n2.vt1)/5;
    time_n4 = sum(n4.vt1)/5;
    time_n8 = sum(n8.vt1)/5;
    
    % Visualization
    f = figure;
    xpoint = [5 10 15 20]; % x-axis
    hold on;
    plot(xpoint,time_n2, '--', "DisplayName","N_{2,0}", "MarkerSize",15);
    plot(xpoint,time_n4, '--', "DisplayName","N_{4,4}", "MarkerSize",15);
    plot(xpoint,time_n8, '--', "DisplayName","N_{8,0}", "MarkerSize",15);
    % change from linear to log in the Y-axis
    set(gca, 'YScale', 'log');
    % Add legend and labels to axis
    legend("Location","best", "FontSize", 15);
    xticks([5 10 15 20]);
    yticks([0.01 0.1 1 10 100]);
    xticklabels([5 10 15 20]);
    yticklabels([0.01 0.1 1 10 100]);
    xlabel('T steps');
    ylabel('Computation Time (s)');
    % Increase font for the paper
    ax = gca;
    ax.XAxis.FontSize = 15; % Set font size of axis
    ax.YAxis.FontSize = 15;
    % save plot
    if is_codeocean
        exportgraphics(f,'/results/logs/rnn_verification_time.pdf', 'ContentType', 'vector');
        saveas(f,'/results/logs/rnn_verification_time.png');
    else
        exportgraphics(f, "../rnn_verification_time.pdf",'ContentType', 'vector');
    end

end

