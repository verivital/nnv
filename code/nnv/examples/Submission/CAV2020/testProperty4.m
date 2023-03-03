function testProperty4()
    % Compare results from CAV 2020 to 2023
    N2 = 1:1:9;
    N1 = 1:1:5;
    results = [];
    results_exact = [];
    
    for n1 = N1
        for n2 = N2
            results = [results; verify_P0_N00_star_appr(n1,n2)];
            results_exact = [results_exact; verify_P0_N00_star(n1, n2)];
        end
    end
    
    save("tP4results", "results", "results_exact");
end