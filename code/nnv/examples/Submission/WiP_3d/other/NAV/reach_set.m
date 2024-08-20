function rT = reach_set()

    %% Reachability analysis of NAV Benchmark
    
    %% Load Components 

    % Load the controller
    % net = importNetworkFromONNX('networks/nn-nav-set.onnx', "InputDataFormats", "BC");
    netonnx = importONNXNetwork('networks/nn-nav-set.onnx', "InputDataFormats", "BC");
    % Load plant
    reachStep = 0.02;
    controlPeriod = 0.2;
    % plant = NonLinearODE(4, 2, @dynamics, reachStep, controlPeriod, eye(4));
    % plant.set_tensorOrder(2);
    % plant.set_taylorTerms(3);
    % plant.set_zonotopeOrder(100);
    % plant.set_intermediateOrder(50);


    %% Reachability analysis

    % Initial set
    lb = [2.9; 2.9; 0; 0];
    ub = [3.1; 3.1; 0; 0];
    init_set = Box(lb,ub);
    init = init_set.partition([1 2],[50 50]);

    % Reachability options
    num_steps = 21;
    reachOptions.reachMethod = 'approx-star';

    N = length(init);
    disp("Verifying "+string(N)+" samples...")
    
    mkdir('temp');
    parpool("Processes"); % initialize parallel process
    
    % Execute reachabilty analysis
    t = tic;
    parfor j = 1:length(init)
        % Get NNV network
        net = matlab2nnv(netonnx);
        % Create plant
        plant = NonLinearODE(4, 2, @dynamics, reachStep, controlPeriod, eye(4));
        plant.set_tensorOrder(2);
        plant.set_taylorTerms(3);
        plant.set_zonotopeOrder(100);
        plant.set_intermediateOrder(50);
        % Get initial conditions
        init_set = init(j).toStar;
        %reachSub = init_set;
        for i = 1:num_steps
            % Compute controller output set
            input_set = net.reach(init_set,reachOptions);

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
    
    f = figure;
    rectangle('Position',[-0.5,-0.5,1,1],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth', 0.1); % goal region
    hold on;
    rectangle('Position',[1,1,1,1],'FaceColor',[0.7 0 0 0.8], 'EdgeColor','r', 'LineWidth', 0.1); % obstacle
    grid;
    t = tic;
    for K = 1 : length(setFiles)
        if ~mod(K,50)
            disp("Plotting partition "+string(K)+" ...");
            toc(t)
            pause(0.01); % to ensure it prints
        end
        res = load("temp/"+setFiles(K).name);
        plant = res.plant;
        % plant.get_interval_sets;
        % Star.plotBoxes_2D_noFill(plant.intermediate_reachSet, 1,2,'b');
        for k=1:(length(plant.cora_set))
            plot(plant.cora_set{k}, [1,2], 'b', 'Unify', true);
        end
    end
    hold on;
    xlabel('x1');
    ylabel('x2');

    disp("Finished plotting all reach sets");


    %% Save figure
    if is_codeocean
        saveas(f,'/results/logs/nav_set.png');
        % exportgraphics(f,'/results/logs/nav-set.pdf', 'ContentType', 'vector');
    else
        saveas(f,'nav_set.png');
        % exportgraphics(f,'nav-set.pdf','ContentType', 'vector');
    end

    % Save results
    if is_codeocean
        save('/results/logs/nav_set.mat','rT','-v7.3');
    else
        save('nav_set.mat', 'rT','-v7.3');
    end

end

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
