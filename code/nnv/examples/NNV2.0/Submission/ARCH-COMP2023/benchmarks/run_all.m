function run_all()

    % Reproduce AINNCS ARCH-COMP 2023 results
    % author: Diego Manzanas
    % date submitted: April XX, 2023
    
    % Supress warnings
    warning ('off','all');
    
    % Turn off figure display
    set(0,'DefaultFigureVisible','off');

%% Begin benchmark evaluation

    % ACC 
    cd ACC;
    % verified (~ 15 seconds)
    acc = reach(); 
    cd ..;
    
    % Benchmark-9 (Tora)
    cd Tora_Heterogeneous;
    % verified with input partition (~ 38 mins)
    tora_relutanh = reachTora_reluTanh();
    % verified with input partition (~ 2.7 hours, let's try a different partition to speed it up)
    tora_sigmoid = reachTora_sigmoid(); 
    cd ..;
    % cd Benchmark9-Tora;
    % reach(); % unknown
    tora_relu = '-';
    % cd ..;
    
    % Benchmark-10 (Unicycle)
    % unknown -> overapproximation
    uncycle = '-';
    
    % VCAS
    cd VCAS;
    % verified (~20 seconds)
    vcas = run_vcas(); 
    cd ..;
    
    % Single Pendulum 
    cd Single_Pendulum;
    % falsified (~2 seconds)
    sp = reach(); 
    cd ..;
    
    % Double Pendulum 
    cd Double_Pendulum;
    run reach_more.m; % falsified
    dp_more = '-';
    run reach_less.m; % falsified
    dp_less = '-';
    cd ..;
    
    % Airplane
    cd Airplane;
    run reach.m; % Falsified
    airplane = '-';
    cd ..;
    
    % Attitude Control
    % Unknown 
    attitude = '-';
    
    % Quad
    % Unknown
    quad = '-';
    
    % Spacecraft Docking 
    % unknown -> overapproximation
    spacecraft = '-';

%% Save results in CSV file for submission

    % Create results variable
    resultsCSV = cell(21, 4); % column names + # of benchmarks (including instances)
    resultsCSV(1, :) = {'benchmark','instance','result','time'}; % columns names

    % ACC
    acc_csv = {'ACC', 'relu', 'verified', acc};
    resultsCSV(2,:) = acc_csv;
    
    % Benchmark-9 (Tora)
    tora_relu_csv = {'TORA', 'relu', '', tora_relu};
    resultsCSV(3,:) = tora_relu_csv;
    tora_relutanh_csv = {'TORA', 'relutanh', '', tora_relutanh};
    resultsCSV(4,:) = tora_relu_csv;
    tora_sigmoid_csv = {'TORA', 'sigmoid', '', tora_sigmoid};
    resultsCSV(5,:) = tora_relu_csv;
    
    % Benchmark-10 (Unicycle)
    
    % VCAS
    
    % Single Pendulum 
    
    % Double Pendulum 
    
    % Airplane
    
    % Attitude Control
    
    % Quad
    
    % Spacecraft Docking 


    % Create and save csv files with results
    if is_codeocean
        writecell(resultsCSV, '/results/logs/results.csv');
    else
        writecell(resultsCSV, 'results.csv');
    end

end