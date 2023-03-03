function HA = buck_v2(~)


%% Generated on 13-May-2020

%-------------Automaton created from Component 'buckboost'-----------------

%% Interface Specification:
%   This section clarifies the meaning of state & input dimensions
%   by showing their mapping to SpaceEx variable names. 

% Component 1 (buckboost.buckboost_template_1):
%  state x := [il; t; gt; D; vc]
%  input u := [uDummy]

%---------------Component buckboost.buckboost_template_1-------------------

%----------------------------State charging--------------------------------

%% equation:
%   il' == (a00c * il + a01c * vc + b0c * Vs) & vc' == (a10c * il + a11c * vc + b1c * Vs) & t' == 1 & gt' == 1 & D' == 0
dynA = ...
[-1288.4615384615381117328070104122,0,0,0,...
-320.51282051282049678775365464389;0,0,0,0,0;0,0,0,0,0;0,0,0,0,0;...
456.62100456621010380331426858902,0,0,0,-45.662100000000002353317540837452];
dynB = ...
[0;0;0;0;0];
dync = ...
[32051.28205128205081564374268055;1;1;0;0];
dynamics = linearSys(dynA, dynB, dync);

%% equation (invariant)
%   t >= 0 & t <= D*T & gt >= 0 & D >= 0 & D <= 1
invA = ...
[0,-1,0,0,0;0,1,0,-1.6667000000000001529599363836454E-05,0;0,0,-1,0,0;0,...
0,0,-1,0;0,0,0,1,0];
invb = ...
[-0;-0;-0;-0;1];

invOpt = struct('A', invA, 'b', invb);
inv = mptPolytope(invOpt);

trans = {};
%% equation (reset)
%   t' == 0
resetA = ...
[1,0,0,0,0;0,0,0,0,0;0,0,1,0,0;0,0,0,1,0;0,0,0,0,1];
resetb = ...
[0;0;0;0;0];
reset = struct('A', resetA, 'b', resetb);

%% equation:
%   t >= D*T
guardA = ...
[0,-1,0,1.6667000000000001529599363836454E-05,0];
guardb = ...
[-0];
guardOpt = struct('A', guardA, 'b', guardb);
guard = mptPolytope(guardOpt);

trans{1} = transition(guard, reset, 2);

loc{1} = location('S1',inv, trans, dynamics);



%---------------------------State discharging------------------------------

%% equation:
%   il' == (a00o * il + a01o * vc + b0o * Vs) & vc' == (a10o * il + a11o * vc + b1o * Vs) & t' == 1 & gt' == 1 & D' == 0
dynA = ...
[-166.66666666666671403618238400668,0,0,0,...
-320.51282051282049678775365464389;0,0,0,0,0;0,0,0,0,0;0,0,0,0,0;...
456.62100456621010380331426858902,0,0,0,-45.662100000000002353317540837452];
dynB = ...
[0;0;0;0;0];
dync = ...
[0;1;1;0;0];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   t >= 0 & t <= (1-D)*T & gt >= 0 & D >= 0 & D <= 1
invA = ...
[0,-1,0,0,0;0,1,0,1.6667000000000001529599363836454E-05,0;0,0,-1,0,0;0,...
0,0,-1,0;0,0,0,1,0];
invb = ...
[-0;1.6667000000000001529599363836454E-05;-0;-0;1];
invOpt = struct('A', invA, 'b', invb);
inv = mptPolytope(invOpt);

trans = {};
%% equation:
%   t' == 0
resetA = ...
[1,0,0,0,0;0,0,0,0,0;0,0,1,0,0;0,0,0,1,0;0,0,0,0,1];
resetb = ...
[0;0;0;0;0];
reset = struct('A', resetA, 'b', resetb);

%% equation:
%   t >= (1-D)*T
guardA = ...
[0,-1,0,-1.6667000000000001529599363836454E-05,0];
guardb = ...
[-1.6667000000000001529599363836454E-05];
guardOpt = struct('A', guardA, 'b', guardb);
guard = mptPolytope(guardOpt);

trans{1} = transition(guard, reset, 1);

loc{2} = location('S2',inv, trans, dynamics);



HA = hybridAutomaton(loc);


end