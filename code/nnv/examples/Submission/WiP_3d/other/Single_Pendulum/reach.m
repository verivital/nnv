% function rT = reach()

    %% Reachability analysis of Single Pendulum Benchmark
    
    %% Load Components 

    % Load the controller
    % net = load_NN_from_mat('controller_single_pendulum.mat');
    % Load plant
    reachStep = 0.01;
    controlPeriod = 0.05;
    % plant = NonLinearODE(3,1,@dynamics_sp, reachStep, controlPeriod, eye(3));
    % plant.set_tensorOrder(2);
    

    %% Reachability analysis

    % Initial set 
    lb = [1.0; 0.0; 0];
    ub = [1.175; 0.2; 0]; % original bounds
    init_set = Box(lb,ub);
    init = init_set.partition([1 2],[35 40]);
    % init_set = Star(lb,ub);
    
    num_steps = 13;
    reachOptions.reachMethod = 'approx-star';

    N = length(init);
    disp("Verifying "+string(N)+" samples...")

    mkdir('temp');
    parpool("Processes"); % initialize parallel process
    
    % Execute reachabilty analysis
    t = tic;
    parfor j = 1:length(init)
        % Create NNCS for each process 
        net = load_NN_from_mat('controller_single_pendulum.mat');
        plant = NonLinearODE(3,1,@dynamics_sp, reachStep, controlPeriod, eye(3));
        plant.set_tensorOrder(2);
        % Get initial conditions
        init_set = init(j).toStar;
        for i = 1:num_steps
            % Compute controller output set
            init_set_s = init_set.affineMap([1 0 0;0 1 0],[]);
            input_set = net.reach(init_set_s,reachOptions);

            % Compute plant reachable set
            init_set = plantReach(plant, init_set, input_set,'lin');
        end
        toc(t);
        parsave("temp/reachSet"+string(j)+".mat",plant);
    end
    rT = toc(t); % get reach time
    disp("Finished reachability...")

    % Shut Down Current Parallel Pool
    poolobj = gcp('nocreate');
    delete(poolobj);

    %% Visualize results
    setFiles = dir('temp/*.mat');
    
    t = tic;
    f = figure;
    hold on;
    rectangle('Position',[0.5,1,1,1],'FaceColor',[1 0 0 0.5],'EdgeColor','r', 'LineWidth',0.1)
    grid;
    for K = 1 : length(setFiles)
        if ~mod(K,50)
            disp("Plotting partition "+string(K)+" ...");
        end
        res = load("temp/"+setFiles(K).name);
        plant = res.plant;
        for k=1:length(plant.cora_set)
            plot(plant.cora_set{k}, [3,1], 'b', 'Unify', true); % this seems to be faster
        end
    end
    xlim([0.1 0.6])
    ylim([0.85 1.25])
    xlabel('Time (s)');
    ylabel('\theta');
    toc(t);

    disp("Finished plotting all reach sets");
    
    %% Save figure
    if is_codeocean
        % exportgraphics(f,'/results/logs/singlePendulum.pdf', 'ContentType', 'vector'); % takes too long
        saveas(f,'/results/logs/singlePendulum.png');
    else
        saveas(f,'singlePendulum.png');
        % exportgraphics(f,'singlePendulum.pdf','ContentType', 'vector');
    end

    % Save results
    if is_codeocean
        save('/results/logs/single_pendulum.mat', 'rT','-v7.3');
    else
        save('single_pendulum.mat', 'rT', '-v7.3');
    end
   
% end

%% Helper function
function init_set = plantReach(plant,init_set,input_set,algoC)
    nS = length(init_set); % based on approx-star, number of sets should be equal
    ss = [];
    for k=1:nS
            ss =[ss plant.stepReachStar(init_set(k), input_set(k),algoC)];
    end
    init_set = ss;
end

function parsave(fname, plant) % trick to save while on parpool
  save(fname, 'plant')
end
