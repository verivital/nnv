function t = reach()

    %% Reachability analysis of Single Pendulum Benchmark
    
    %% Load Components 

    % Load the controller
    net = load_NN_from_mat('controller_single_pendulum.mat');
    % Load plant
    reachStep = 0.01;
    controlPeriod = 0.05;
    plant = NonLinearODE(3,1,@dynamics_sp, reachStep, controlPeriod, eye(3));
    
    
    %% Random simulations (for falsification)

    rng(0);
    % Initial region
    lb = [1.0; 0.0]; 
    ub = [1.2; 0.2];
    N = 1000; % number of simulations
    num_steps = 11;% initial + 10 control steps
    r1 = unifrnd(lb(1),ub(1),[N 1]);
    r2 = unifrnd(lb(2),ub(2),[N 1]);
    r = [r1 r2 zeros(N,1)];
    sims = zeros(3,N,num_steps);
    t = tic;
    for k=1:N
        rK = r(k,:); % initial state
        sims(:,k,1) = rK(end,:);
        cA = 0; % control action
        for s = 2:num_steps+1
            cA = net.evaluate(rK(end,1:2)');
            [~,rK] = plant.evaluate(rK(end,:), cA);
            sims(:,k,s) = rK(end,:);
        end
        if sims(1,k,11) >= 1
            trace = k;
            break;
        end
    end
    t = toc(t);
    
    %% Visualize results
    
    f = figure;
    hold on;
    rectangle('Position',[0.5,1,1,1],'FaceColor',[1 0 0 0.5],'EdgeColor','r', 'LineWidth',0.1)
    % Plot only falsifying trace
    plot(squeeze(sims(3,k,:)), squeeze(sims(1,k,:)), 'Color', [0 0 1 1]);
    grid;
    xlabel('Time (s)');
    ylabel('\theta');
    xlim([0 0.6])
    ylim([0.95 1.25])
    % Save figure
    if is_codeocean
        exportgraphics(f,'/results/logs/singlePendulum.pdf', 'ContentType', 'vector');
    else
        exportgraphics(f,'singlePendulum.pdf','ContentType', 'vector');
    end

end