function run_all()
    
    % Supress warnings
    warning ('off','all');
    
    % Turn off figure display
    set(0,'DefaultFigureVisible','off');


%% Begin benchmark evaluation

    % ACC 
    
    cd ACC;
    
    try
        acc = reach(); % verified (~ 15 seconds)
    catch
        acc = -1;
    end
    
    cd ..;

    % Benchmark-9 (Tora)
    cd Tora_Heterogeneous;
    
    try
        tora_relutanh = reachTora_reluTanh(); % verified w/ input partition (~ 38 mins)
    catch
        tora_relutanh = -1;
    end
    
    try
        tora_sigmoid = reachTora_sigmoid();   % verified w/ input partition (~ 1 hour)
    catch
        tora_sigmoid = -1;
    end
    
    cd ..;

    cd Benchmark9-Tora;
    
    try
        tora_relu = reach(); % verified (~25 seconds)
    catch
        tora_relu = -1;
    end
    
    cd ..;

    % Benchmark-10 (Unicycle)
    
    cd Benchmark10-Unicycle/;
    
    try
        unicycle = reach(); % overapprox -> unknown
    catch
        unicycle = -1;
    end
    
    % VCAS
    
    cd VCAS;
    
    try
        vcas = run_vcas(); % verified (~20 seconds total)
    catch
        vcas = -1;
    end
    
    cd ..;

    % Single Pendulum 
    
    cd Single_Pendulum;
    
    try 
        sp = reach(); % falsified (reach sets ~3 seconds)
    catch
        sp = -1;
    end
    
    cd ..;

    % Double Pendulum 
    
    cd Double_Pendulum;
    
    try
        dp_more = reach_more(); % falsified (reach sets ~ 30 seconds)
    catch
        dp_more = -1;
    end

    try
        dp_less = reach_less(); % falsified (reach sets ~ 30 seconds)
    catch
        dp_less = -1;
    end
    
    cd ..;

    % Airplane
    
    cd Airplane;
    
    try
        airplane = reach(); % Falsified (reach sets ~ 7 seconds)
    catch
        airplane = -1;
    end
    
    cd ..;
    
    % Attitude Control
    
    cd('Attitude Control');
    
    try
        attitude = reach(); % sigmoid controller -> overapprox -> unknown
    catch 
        attitude = -1;
    end
    
    cd ..;
    
    % Quad
    
    cd QUAD;
    
    try 
        quad = reach();
        % quad = '-'; % sigmoid controller -> overapprox -> unknown
    catch
        quad = -1;
    end
    
    cd ..;
    
    % Spacecraft Docking 
    
    cd Docking;
    
    try
        spacecraft = reach(); % overapproximation -> unknown
    catch
        spacecraft = -1;
    end
    
    cd ..;

    %%%%%%%%% New from 2024 %%%%%%%%%%%

    % NAV
    
    cd NAV/;
    
    try
        navPoint = reach_point();
    catch
        navPoint = -1;
    end
    
    try
        navSet = reach_set();
    catch
        navSet = -1;
    end
    
    % Cartpole
    
    cd Cartpole;
    
    try
        cartpole = reach();
    catch
        cartpole = -1;
    end

%% Save results in CSV file for submission

    % Create results variable
    resultsCSV = cell(24, 4); % column names + # of benchmarks (including instances)
    resultsCSV(1, :) = {'benchmark','instance','result','time'}; % columns names

    % ACC
    acc_csv = {'ACC', 'safety', 'verified', acc};
    resultsCSV(2,:) = acc_csv;

    % Airplane
    airplane_csv = {'Airplane', 'continuous', 'violated', airplane};
    resultsCSV(3,:) = airplane_csv;


    % Attitude Control
    attitude_csv = {'AttitudeControl','avoid','unknown',attitude};
    resultsCSV(4,:) = attitude_csv;

    % Cartpole
    cartpole_csv = {'Cartpole', 'reach', 'unknown', cartpole};
    resultsCSV(5,:) = cartpole_csv;

    % Spacecraft Docking 
    docking_csv = {'Docking', 'constraint', 'unknown', spacecraft};
    resultsCSV(6,:) = docking_csv;

    % Double Pendulum
    dp_csv(1,:) = {'DoublePendulum', 'more robust', 'violated', dp_more};
    dp_csv(2,:) = {'DoublePendulum', 'less robust', 'violated', dp_less};
    resultsCSV(7:8, :) = dp_csv;

    % NAV
    dp_csv(1,:) = {'NAV', 'standard', '?', navPoint};
    dp_csv(2,:) = {'NAV', 'robust', '?', navSet};
    resultsCSV(9:10, :) = dp_csv;

    % Quad
    quad_csv = {'QUAD', 'reach', 'unknown', quad};
    resultsCSV(11,:) = quad_csv;

    % Single Pendulum 
    sp_csv = {'SinglePendulum', 'reach', 'violated', sp};
    resultsCSV(12,:) = sp_csv; 

    % Benchmark-9 (Tora)
    tora_csv(1,:) = {'TORA', 'remain',         'verified', tora_relu};
    tora_csv(2,:) = {'TORA', 'reach-tanh',     'verified', tora_relutanh};
    tora_csv(3,:) = {'TORA', 'reach-sigmoid',  'verified', tora_sigmoid};
    resultsCSV(13:15,:) = tora_csv;

    % Benchmark-10 (Unicycle)
    unicycle_csv = {'Unicycle', 'reach', 'unknown', unicycle};
    resultsCSV(16,:) = unicycle_csv;

    % VCAS
    vcas_m19 = {'VCAS', 'middle19', 'verified', vcas(1,1)};
    vcas_m22 = {'VCAS', 'middle22', 'verified', vcas(2,1)};
    vcas_m25 = {'VCAS', 'middle25', 'violated', vcas(3,1)};
    vcas_m28 = {'VCAS', 'middle28', 'violated', vcas(4,1)};
    vcas_w19 = {'VCAS', 'worst19',  'verified', vcas(1,2)};
    vcas_w22 = {'VCAS', 'worst22',  'violated', vcas(2,2)};
    vcas_w25 = {'VCAS', 'worst25',  'violated', vcas(3,2)};
    vcas_w28 = {'VCAS', 'worst28',  'violated', vcas(4,2)};
    resultsCSV(17:24,:) = [vcas_m19; vcas_m22; vcas_m25; vcas_m28; vcas_w19; vcas_w22; vcas_w25; vcas_w28];


    % Create and save csv files with results
    if is_codeocean
        writecell(resultsCSV, '/results/results.csv'); % ARCH repeatability specific
    else
        writecell(resultsCSV, 'results.csv');
    end

end