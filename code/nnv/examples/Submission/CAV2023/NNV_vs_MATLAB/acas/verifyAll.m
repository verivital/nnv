%% Verify properties 1 to 4 of ACAS Xu (All networks)

%% Define reachability options
reachOptionsList = [];
% approx star
reachOpt = struct;
reachOpt.reachMethod = 'approx-star';
reachOpt.relaxFactor = 0;
reachOptionsList = [reachOptionsList; reachOpt];
% relax star
reachOpt.relaxFactor = 0.2;
reachOptionsList = [reachOptionsList; reachOpt];
reachOpt.relaxFactor = 0.4;
reachOptionsList = [reachOptionsList; reachOpt];
reachOpt.relaxFactor = 0.6;
reachOptionsList = [reachOptionsList; reachOpt];
reachOpt.relaxFactor = 0.8;
reachOptionsList = [reachOptionsList; reachOpt];
reachOpt.relaxFactor = 1;
reachOptionsList = [reachOptionsList; reachOpt];
% exact-star
% reachOpt = struct;
% reachOpt.reachMethod = 'exact-star';
% reachOpt.relaxFactor = 0; % dummy variable, just to match dimensions of other reachOpt, not used in reachability
% reachOptionsList = [reachOptionsList; reachOpt];

%% Verify Properties

% Property 1
verifyP1(reachOptionsList);

% Property 2
verifyP2(reachOptionsList);

% Property 3
verifyP3(reachOptionsList);

% Property 4
verifyP4(reachOptionsList);
