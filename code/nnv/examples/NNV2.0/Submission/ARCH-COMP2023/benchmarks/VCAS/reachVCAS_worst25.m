function t = reachVCAS_worst25()    
    %% Reachability analysis of VCAS
    
    % 1) Load components and set reachability parameters

    % Load controllers
    networks.vcas1 = load_NN_from_mat('nnv_networks/VertCAS_noResp_pra01_v9_20HU_200.mat');
    networks.vcas2 = load_NN_from_mat('nnv_networks/VertCAS_noResp_pra02_v9_20HU_200.mat');
    networks.vcas3 = load_NN_from_mat('nnv_networks/VertCAS_noResp_pra03_v9_20HU_200.mat');
    networks.vcas4 = load_NN_from_mat('nnv_networks/VertCAS_noResp_pra04_v9_20HU_200.mat');
    networks.vcas5 = load_NN_from_mat('nnv_networks/VertCAS_noResp_pra05_v9_20HU_200.mat');
    networks.vcas6 = load_NN_from_mat('nnv_networks/VertCAS_noResp_pra06_v9_20HU_200.mat');
    networks.vcas7 = load_NN_from_mat('nnv_networks/VertCAS_noResp_pra07_v9_20HU_200.mat');
    networks.vcas8 = load_NN_from_mat('nnv_networks/VertCAS_noResp_pra08_v9_20HU_200.mat');
    networks.vcas9 = load_NN_from_mat('nnv_networks/VertCAS_noResp_pra09_v9_20HU_200.mat');
    
    % Controller outputs
    % 1) COC
    % 2) DNC
    % 3) DND
    % 4) DES1500
    % 5) CL1500
    % 6) SDES1500
    % 7) SCL1500
    % 8) SDES2500
    % 9) SCL2500
    
    % Accelerations
    g = 32.2; %ft/s^2
    h0.COC = [-g/8, 0, g/8];
    h0.DNC = [-g/3, -7*g/24, -g/4];
    h0.DND = [g/4, 7*g/24, g/3];
    h0.DES1500 = [-g/3, -7*g/24, -g/4];
    h0.CL1500 = h0.DND;
    h0.SDES1500 = -g/3;
    h0.SCL1500 = g/3;
    h0.SDES2500 = -g/3;
    h0.SCL2500 = g/3;
    
    % Create plant
    controlPeriod = 1;
    out_mat = [1 0 0 0;0 1 0 0;0 0 1 0];
    plant = DNonLinearODE(4,2,@planeDynamics, controlPeriod, out_mat);
    
    %% Reachability analysis

    % Set parameters and options for reachability
    
    % Initial set
    lb = [-133; -25.5; 25; 1];
    ub = [-129; -25.5; 25; 1];
    init_set = Star(lb,ub); % initial state set
    minIdx = 1;
    reachOptions.reachMethod = 'approx-star';

    % initialize variables
    reachAll = init_set;
    steps = 10;
    uNN_all = cell(1,steps);
    yNN_all = cell(1,steps);
    idx_all = cell(1,steps);
    inp_all = cell(1,steps);
    rand_idx = cell(1,steps);

    % Execute reachabilty analysis
    t = tic;
    for i=1:steps
        % Compute plant output set given states
        uNN = plantOut(init_set,out_mat);
        % Compute controller output set
        yNN = reach_nn(minIdx, uNN, networks, reachOptions);
        % Get advisory value
        minIdx = getMaxIndexes(yNN);
        [rand_idx{i},input_set] = advisoryVCAS(minIdx,h0,init_set);
        % Store sets and advisories at each control step
        uNN_all{i} = uNN;
        yNN_all{i} = yNN;
        idx_all{i} = minIdx;
        inp_all{i} = input_set;
        % Compute plant reachable set
        init_set = plantReach(plant, init_set, input_set);
        % store state trajectory
        reachAll = [reachAll init_set]; 
    end
    t = toc(t);

    % Save results
    if is_codeocean
        save('/results/logs/VCAS_worst25.mat', 'reachAll','t','-v7.3');
    else
        save('VCAS_worst25.mat', 'reachAll','t','-v7.3');
    end


    %% Visualize results

    times = 0:controlPeriod:steps*controlPeriod;
    f = figure;
    rectangle('Position',[0,-100,10,200],'FaceColor',[0.5 0 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
    hold on;
    Star.plotRanges_2D(reachAll,1,times,'b');
    grid;
    xlabel('Time (s)');
    ylabel('Distance (ft)');
    
    % Save figure
    if is_codeocean
        exportgraphics(f,'/results/logs/VCAS_worst25.pdf', 'ContentType', 'vector');
    else
        exportgraphics(f,'VCAS_worst25.pdf','ContentType', 'vector');
    end


end



%% Helper Functions

% Choose NN to execute
function y = reach_nn(advise,Unn,networks,method)
    y = [];
    scale_mean = [0.0,0.0,20.0,-0.432599379632]; % Means    
    scale_range =[16000.0,5000.0,40.0,3.102300001]; % Ranges
    Unn = Unn.affineMap(eye(3),-1*scale_mean(1:3)');
    Unn = Unn.affineMap(diag(1./scale_range(1:3)),[]);
    for a=1:length(advise)
        if advise(a) == 1
            y = [y networks.vcas1.reach(Unn,method)];
        elseif advise(a) == 2
            y = [y networks.vcas2.reach(Unn,method)];
        elseif advise(a) == 3
            y = [y networks.vcas3.reach(Unn,method)];
        elseif advise(a) == 4
            y = [y networks.vcas4.reach(Unn,method)];
        elseif advise(a) == 5
            y = [y networks.vcas5.reach(Unn,method)];
        elseif advise(a) == 6
            y = [y networks.vcas6.reach(Unn,method)];
        elseif advise(a) == 7
            y = [y networks.vcas7.reach(Unn,method)];
        elseif advise(a) == 8
            y = [y networks.vcas8.reach(Unn,method)];
        elseif advise(a) == 9
            y = [y networks.vcas9.reach(Unn,method)];
        end
    end
end

% Set advisory to ownship
function [idx,y] = advisoryVCAS(r,accs, uNN)
    y = [];
    idx = 1;
    [m,M] = uNN.getRanges;
    for l=1:length(r)
        if r(l) == 1
            y = [y Star([1;accs.COC(idx)],[1;accs.COC(idx)])];
        elseif r(l) == 2
            if M(2) < 0 && (M(4) == 2 || m(4) == 2)
                y = [y Star([2;0],[2;0])];
            else
                y = [y Star([2;accs.DNC(idx)],[2;accs.DNC(idx)])];
            end
        elseif r(l) == 3
            if m(2) > 0 && (M(4) == 3 || m(4) == 3)
                y = [y Star([3;0],[3;0])];
            else
                y = [y Star([3;accs.DND(idx)],[3;accs.DND(idx)])];
            end
        elseif r(l) == 4
            if M(2) < -25 && (M(4) == 4 || m(4) == 4)
                y = [y Star([4;0],[4;0])];
            else
                y = [y Star([4;accs.DES1500(idx)],[4;accs.DES1500(idx)])];
            end
        elseif r(l) == 5
            if m(2) > 25 && (M(4) == 5 || m(4) == 5)
                y = [y Star([5;0],[5;0])];
            else
                y = [y Star([5;accs.CL1500(idx)],[5;accs.CL1500(idx)])];
            end
        elseif r(l) == 6
            if M(2) < -25&& (M(4) == 6 || m(4) == 6)
                y = [y Star([6;0],[6;0])];
            else
                y = [y Star([6;accs.SDES1500],[6;accs.SDES1500])];
            end
        elseif r(l) == 7
            if m(2) > 25 && (M(4) == 7 || m(4) == 7)
                y = [y Star([7;0],[7;0])];
            else
                y = [y Star([7;accs.SCL1500],[7;accs.SCL1500])];
            end
        elseif r(l) == 8
            if M(2) < -41.66 && (M(4) == 8 || m(4) == 8)
                y = [y Star([8;0],[8;0])];
            else
                y = [y Star([8;accs.SDES2500],[8;accs.SDES2500])];
            end
        elseif r(l) == 9
            if m(2) > 41.66 && (M(4) == 9 || m(4) == 9)
                y = [y Star([9;0],[9;0])];
            else
                y = [y Star([9;accs.SCL2500],[9;accs.SCL2500])];
            end
        end
    end
end

% Get min index for AcasXu networks (advisory command)
function idxs = getMaxIndexes(star_set)
    if length(star_set) > 1
        Xi = Star.merge_stars(star_set,1,'single');
    else
        Xi = star_set;
    end
    [mm,MM] = Xi.getRanges;
    [~, idx1] = max(mm);
    idxs = [];
    for idxm=1:length(mm)
        if MM(idxm) >= mm(idx1)
            idxs = [idxs idxm];
        end
    end
end

% Compute output set of plant
function Rout = plantOut(init_set,out_mat)
    Rout = [];
    for o=1:length(init_set)
        Rout =[Rout init_set(o).affineMap(out_mat,[])];
    end
end

% Compute stat sets for the plant
function stateSet = plantReach(plant, init_set, Up)
    stateSet = [];
    for s=1:length(init_set)
        for u=1:length(Up)
            stateSet = [stateSet plant.stepReachStar(init_set(s),Up(u))];
        end
    end    
end