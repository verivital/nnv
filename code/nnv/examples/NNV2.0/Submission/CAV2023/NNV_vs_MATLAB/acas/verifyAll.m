%% Verify properties 1 to 4 of ACAS Xu (All networks)
function verifyAll()

    %% Define reachability options

    reachOptionsList = [];

    % approx star
    reachOpt = struct;
    reachOpt.reachOption = "single";
    reachOpt.numCores = 1;
    reachOpt.reachMethod = 'approx-star';
    reachOpt.relaxFactor = 0;
    reachOptionsList = [reachOptionsList; reachOpt]; % approx

    % relax star
    reachOpt.reachMethod = 'relax-star-range';
    reachOpt.relaxFactor = 0.25;
    reachOptionsList = [reachOptionsList; reachOpt]; % relax 0.25
    reachOpt.relaxFactor = 0.5;
    reachOptionsList = [reachOptionsList; reachOpt]; % relax 0.5
    reachOpt.relaxFactor = 0.75;
    reachOptionsList = [reachOptionsList; reachOpt]; % relax 0.75
    reachOpt.relaxFactor = 1;
    reachOptionsList = [reachOptionsList; reachOpt]; % relax 1

    % exact-star
    reachOpt.reachMethod = 'exact-star';
    reachOpt.reachOption = "parallel";
    max_cores = maxNumCompThreads;
    reachOpt.numCores = min([8,max_cores]);
    reachOptionsList = [reachOptionsList; reachOpt]; % exact (parallel)
    
    %% Verify Properties
    
    % Property 3
    verifyP3(reachOptionsList);

    % Property 4
    verifyP4(reachOptionsList);

    % Close parallel pool if still active
    poolobj = gcp('nocreate');
    delete(poolobj);

end