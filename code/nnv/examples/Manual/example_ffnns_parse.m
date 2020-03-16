% /* Example of parsing a network trained from Matlab */
load Engine_Toy_Tansig_net.mat; % load the network 
F = FFNNS.parse(net); % parse the network
