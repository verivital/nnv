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
    reachOptionsList = [reachOptionsList; reachOpt]; % relax
    reachOpt.relaxFactor = 0.5;
    reachOptionsList = [reachOptionsList; reachOpt]; % relax
    reachOpt.relaxFactor = 0.75;
    reachOptionsList = [reachOptionsList; reachOpt]; % relax
    reachOpt.relaxFactor = 1;
    reachOptionsList = [reachOptionsList; reachOpt]; % relax
    % exact-star
    reachOpt.reachMethod = 'exact-star';
%     reachOptionsList = [reachOptionsList; reachOpt]; % exact (single)
    % exact-star (parallel)
    reachOpt.reachOption = "parallel";
    max_cores = getenv('NUMBER_OF_PROCESSORS');
    reachOpt.numCores = min(8,max_cores);
    reachOptionsList = [reachOptionsList; reachOpt]; % exact (parallel)
    
    %% Verify Properties
    
    % Property 3
    verifyP3(reachOptionsList);
    % 
    % % Property 4
    verifyP4(reachOptionsList);
end