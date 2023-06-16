function run_all()

    % Reproduce AINNCS ARCH-COMP 2023 results
    % author: Diego Manzanas
    % date submitted: April 19, 2023
    % 
    % TODO: prepare dockerfile for repeatability
    
    % Supress warnings
    warning ('off','all');
    
    % Turn off figure display
    set(0,'DefaultFigureVisible','off');

%% Begin benchmark evaluation

    % ACC 
    cd ACC;
    acc = reach(); % verified (~ 15 seconds)
    cd ..;
    
    % Benchmark-9 (Tora)
    cd Tora_Heterogeneous;
    tora_relutanh = reachTora_reluTanh(); % verified w/ input partition (~ 38 mins)
    tora_sigmoid = reachTora_sigmoid();   % verified w/ input partition (~ 1 hour)
    cd ..;
    cd Benchmark9-Tora;
    tora_relu = reach(); % verified (~25 seconds)
    cd ..;
    
    % Benchmark-10 (Unicycle)
    unicycle = '-'; % overapprox -> unknown
    
    % VCAS
    cd VCAS;
    vcas = run_vcas(); % verified (~20 seconds)
    cd ..;
    
    % Single Pendulum 
    cd Single_Pendulum;
    sp = reach(); % falsified (reach sets ~3 seconds)
    cd ..;
    
    % Double Pendulum 
    cd Double_Pendulum;
    dp_more = reach_more(); % falsified (reach sets ~ 30 seconds)
    dp_less = reach_less(); % falsified (reach sets ~ 30 seconds)
    cd ..;
    
    % Airplane
    cd Airplane;
    airplane = reach(); % Falsified (reach sets ~ 7 seconds)
    cd ..;
    
    % Attitude Control
    attitude = '-'; % sigmoid controller -> overapprox -> unknown
    
    % Quad
    quad = '-'; % sigmoid controller -> overapprox -> unknown
    
    % Spacecraft Docking 
    spacecraft = '-'; % overapproximation -> unknown

%% Save results in CSV file for submission

    % Create results variable
    resultsCSV = cell(21, 4); % column names + # of benchmarks (including instances)
    resultsCSV(1, :) = {'benchmark','instance','result','time'}; % columns names

    % ACC
    acc_csv = {'ACC', 'relu', 'verified', acc};
    resultsCSV(2,:) = acc_csv;
    
    % Benchmark-9 (Tora)
    tora_csv(1,:) = {'TORA', 'relu',     'verified', tora_relu};
    tora_csv(2,:) = {'TORA', 'relutanh', 'verified', tora_relutanh};
    tora_csv(3,:) = {'TORA', 'sigmoid',  'verified', tora_sigmoid};
    resultsCSV(3:5,:) = tora_csv;
    
    % Benchmark-10 (Unicycle)
    unicycle_csv = {'Unicycle', '', 'unknown', unicycle};
    resultsCSV(6,:) = unicycle_csv;
    
    % VCAS
    vcas_m19 = {'VCAS', 'middle19', 'verified', vcas(1,1)};
    vcas_m22 = {'VCAS', 'middle22', 'verified', vcas(2,1)};
    vcas_m25 = {'VCAS', 'middle25', 'violated', vcas(3,1)};
    vcas_m28 = {'VCAS', 'middle28', 'violated', vcas(4,1)};
    vcas_w19 = {'VCAS', 'worst19',  'verified', vcas(1,2)};
    vcas_w22 = {'VCAS', 'worst22',  'violated', vcas(2,2)};
    vcas_w25 = {'VCAS', 'worst25',  'violated', vcas(3,2)};
    vcas_w28 = {'VCAS', 'worst28',  'violated', vcas(4,2)};
    resultsCSV(7:14,:) = [vcas_m19; vcas_m22; vcas_m25; vcas_m28; vcas_w19; vcas_w22; vcas_w25; vcas_w28];
    
    % Single Pendulum 
    sp_csv = {'SinglePendulum', '', 'violated', sp};
    resultsCSV(15,:) = sp_csv; 
    
    % Double Pendulum
    dp_csv(1,:) = {'DoublePendulum', 'more', 'violated', dp_more};
    dp_csv(2,:) = {'DoublePendulum', 'less', 'violated', dp_less};
    resultsCSV(16:17, :) = dp_csv;
    
    % Airplane
    airplane_csv = {'Airplane', '', 'violated', airplane};
    resultsCSV(18,:) = airplane_csv;
    
    % Attitude Control
    attitude_csv = {'AttitudeControl','','unknown',attitude};
    resultsCSV(19,:) = attitude_csv;
    
    % Quad
    quad_csv = {'Quadrotor', '', 'unknown', quad};
    resultsCSV(20,:) = quad_csv;
    
    % Spacecraft Docking 
    spacecraft_csv = {'Spacecraft', '','unknown',spacecraft};
    resultsCSV(21,:) = spacecraft_csv;

    
    
    % Create and save csv files with results
    if is_codeocean
        writecell(resultsCSV, '/results/results.csv'); % ARCH repeatability specific
    else
        writecell(resultsCSV, 'results.csv');
    end

end