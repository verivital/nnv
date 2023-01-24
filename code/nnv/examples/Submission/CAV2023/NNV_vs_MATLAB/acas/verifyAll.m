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
    reachOpt.relaxFactor = 0.25;
    reachOptionsList = [reachOptionsList; reachOpt]; % relax
    reachOpt.relaxFactor = 0.5;
    reachOptionsList = [reachOptionsList; reachOpt]; % relax
    reachOpt.relaxFactor = 0.75;
    reachOptionsList = [reachOptionsList; reachOpt]; % relax
    reachOpt.relaxFactor = 1;
    reachOptionsList = [reachOptionsList; reachOpt]; % relax
    reachOpt.reachMethod = 'exact-star';
    reachOptionsList = [reachOptionsList; reachOpt]; % exact (single)
    reachOpt.reachOption = "parallel";
    reachOpt.numCores = 8;
    reachOptionsList = [reachOptionsList; reachOpt]; % exact (parallel)
    
    %% Verify Properties
    
    % Property 1
%     verifyP1(reachOptionsList);
    
    % Property 2
%     verifyP2(reachOptionsList);
    
    % Property 3
%     verifyP3(reachOptionsList);
    % 
    % % Property 4
    verifyP4(reachOptionsList);
end