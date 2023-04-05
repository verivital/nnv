function reach_acc_linear()

    % Reachability analysis of Adaptive Cruise Control with a neuralODE (linear) as a plant model

    % safety specification: relative distance > safe distance
    % dis = x_lead - x_ego  
    % safe distance between two cars, see here 
    % https://www.mathworks.com/help/mpc/examples/design-an-adaptive-cruise-control-system-using-model-predictive-control.html
    % dis_safe = D_default + t_gap * v_ego;

    %%  Load objects
    
    % Load NN
    net = load_NN_from_mat('controller_5_20.mat');

    % Load plant data
    plant_data = load('plant_3rd_order_node.mat');
    states = 8; 
    
    % Extract neural ode data
    w1 = plant_data.neuralOdeParameters.fc1.Weights;
    b1 = plant_data.neuralOdeParameters.fc1.Bias;
    w2 = plant_data.neuralOdeParameters.fc2.Weights;
    b2 = plant_data.neuralOdeParameters.fc2.Bias;
    % Transform neuralODE data into linearODE format
    A = w2*w1;
    B = b2 + w2*b1;
    A3rd = zeros(8,8);
    A3rd(1,:) = [0 1 0 0 0 0 0 0];
    A3rd(2,:) = [0 0 1 0 0 0 0 0];
    A3rd(3,:) = [0 0 A(1,1) 0 0 A(1,2:4)];
    A3rd(4,:) = [0 0 0 0 1 0 0 0];
    A3rd(5,:) = [0 0 0 0 0 1 0 0];
    A3rd(6,:) = [0 0 A(2,1) 0 0 A(2,2:4)];
    A3rd(7,:) = [0 0 A(3,1) 0 0 A(3,2:4)];
    A3rd(8,:) = [0 0 A(4,1) 0 0 A(4,2:4)];
    B3rd = extractdata([zeros(2,1);B(1);zeros(2,1);B(2:4)]);
    C = eye(states); C(7,7) = 0; C(end) = 0;
    D = zeros(8,1);
    cp = 0.1;
    
    % Create plant model from neural ode
    plant = LinearODE(A3rd,B3rd,C,D,cp,10);
    

    %% Reachability analysis
    
    % Set reachability options

    % initial condition of x_lead
    xlead = [90 110];
    % initial condition of v_lead
    v_lead = [32 32.2];
    % initial condition of x_ego
    x_ego = [10 11]; 
    % initial condition of v_ego
    v_ego = [30 30.2];
    % input set
    U = Star(0,0);
    % Initial state set
    % [xlead,vlead,a_int_lead,xego,vego,z_int_ego,aego,alead]
    lb = [xlead(1); v_lead(1); 0; x_ego(1); v_ego(1); 0; 0; -2]; %lower bound
    ub = [xlead(2); v_lead(2); 0; x_ego(2); v_ego(2); 0; 0; -2]; % upper bound
    X0 = Star(lb,ub); %
    map_mat = [0 0 0 0 1 0 0 0;
                1 0 0 -1 0 0 0 0;
                0 1 0 0 -1 0 0 0];
    U_fix = Star([30;1.4],[30;1.4]); % vset and tgap
    
    % Perform reachability analysis
    trajR = X0;
    trajU = [];
    N = 50;
    ts = 0.01;
    % Start reachability
    t = tic;
    for k=1:N
        R1 = plant.simReach('direct', X0, U, ts, 10, []); % reachability of plant
        X0 = R1(end); % Get only last reach set (reach set at control period)
        trajR = [trajR X0]; % Keep track of trajectory of NNCS
        ppp = X0.affineMap(map_mat,[]);
        Uin = U_fix.concatenate(ppp);
        Rc = net.reach(Uin); % controller reach comp
        trajU = [trajU Rc];
        x08 = X0.affineMap([0 0 0 0 0 0 0 1],[]);
        X0 = X0.affineMap(C(1:6,:),[]); % Get set for variables 1 to 6
        X0 = X0.concatenate(Rc); % Add state/input 7 (a_ego)
        X0 = X0.concatenate(x08);
        if X0.nVar > 100
            X0 = X0.getBox;
            X0 = X0.toStar;
        end
    end
    rT = toc(t);

    % Save results
    save('results_linearplant.mat','rT',"trajR","trajU");
    

    %% Visualize results

    % Transform reach sets for visualization
    t_gap = 1.4;
    D_default = 10;
    alp = 1;
    outAll = [];
    safe_dis = [];
    for i=1:length(trajR)
        outAll = [outAll trajR(i).affineMap([1 0 0 -1 0 0 0 0], [])]; % distance betweeen vehicles
        safe_dis = [safe_dis trajR(i).affineMap([0 0 0 0 alp*t_gap 0 0 0], alp*D_default)]; % safety distance
    end
    times = 0:0.1:0.1*N; % to plot in x-axis (time)

    % Create figure
    f = figure;
    hold on;
    pb = plot(0,85,'b');
    pr = plot(0,85,'m');
    Star.plotRanges_2D(outAll,1,times,'b');
    hold on;
    Star.plotRanges_2D(safe_dis,1,times,'m');

    % Enhance figure for paper
    ax = gca; % Get current axis
    ax.XAxis.FontSize = 15; % Set font size of axis
    ax.YAxis.FontSize = 15;
    xlabel('Time (s)');
    ylabel('Distance (m)')
    legend([pr,pb],{'safe dist','rel dist'},"Location","best",'FontSize',14);

    % Save results
    if is_codeocean
        exportgraphics(f,'/results/logs/acc_linear.pdf','ContentType','vector');
    else
        exportgraphics(f,'acc_linear.pdf','ContentType','vector');
    end

end