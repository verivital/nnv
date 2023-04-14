function time = run_vcas()
    time = zeros(4,2);
    % VCAS middle acceleration
    time(1,1) = reachVCAS_middle19();
    time(2,1) = reachVCAS_middle22();
    time(3,1) = reachVCAS_middle25();
    time(4,1) = reachVCAS_middle28();
    % VCAS worst acceleration
    time(1,2) = reachVCAS_worst19();
    time(2,2) = reachVCAS_worst22();
    time(3,2) = reachVCAS_worst25();
    time(4,2) = reachVCAS_worst28();
end