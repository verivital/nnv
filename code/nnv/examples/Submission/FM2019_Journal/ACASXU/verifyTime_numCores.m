numCores = [1 2 3 4 5 6];
n1 = 1;
n2 = 3; 
property_id = 4;
method = 'exact-star';
N = length(numCores);
vt = zeros(1, N);
for i=1:N
    [~,vt(i)] = verifyNetwork(n1, n2, property_id, method, numCores(i));
end

% plot verification time versus numcores
ReluplexTime = 1150*ones(1, N);
fig = figure;
plot(numCores, vt, '-*', 'LineWidth',3);
hold on;
plot(numCores, ReluplexTime, '-x', 'LineWidth',3);
xlabel('N', 'Fontsize', 16);
ylabel('VT (sec)', 'Fontsize', 16);
legend('Star', 'Reluplex');
set(gca, 'Fontsize', 16);
saveas(gcf, 'verificationTime_vs_numCores.pdf');


%% function
function [safe, vt] = verifyNetwork(n1, n2, property_id, method, numCores)
    [F, I] = constructNetwork_inputSet(n1, n2, property_id);
    U = getSpecs(property_id);
    [safe, vt, ~] = F.verify(I, U, method, numCores);
end


% construct a network from network id
function [F, I] = constructNetwork_inputSet(n1,n2,property_id)    
    addpath(genpath("../../../examples/NN/ACASXU/nnet-mat-files/"));
    load(['ACASXU_run2a_',num2str(n1),'_',num2str(n2),'_batch_2000.mat']); 
    Layers = [];
    n = length(b);
    for i=1:n - 1
        bi = cell2mat(b(i));
        Wi = cell2mat(W(i));
        Li = LayerS(Wi, bi, 'poslin');
        Layers = [Layers Li];
    end
    bn = cell2mat(b(n));
    Wn = cell2mat(W(n));
    Ln = LayerS(Wn, bn, 'purelin');
    Layers = [Layers Ln];
    F = FFNNS(Layers);
    
    [lb, ub] = getInputs(property_id);
    % normalize input
    for i=1:5
        lb(i) = (lb(i) - means_for_scaling(i))/range_for_scaling(i);
        ub(i) = (ub(i) - means_for_scaling(i))/range_for_scaling(i);   
    end
    
    I = Star(lb, ub);
        
end

function [lb, ub] = getInputs(property_id)
    if property_id == 1 || property_id == 2
        % Input Constraints
        % 55947.69 <= i1(\rho) <= 60760,
        % -3.14 <= i2 (\theta) <= 3.14,
        %-3.14 <= i3 (\shi) <= -3.14
        % 1145 <= i4 (\v_own) <= 1200, 
        % 0 <= i5 (\v_in) <= 60
        lb = [55947.69; -3.14; -3.14; 1145; 0];
        ub = [60760; 3.14; 3.14; 1200; 60];

    elseif property_id == 3
        % Input Constraints
        % 1500 <= i1(\rho) <= 1800,
        % -0.06 <= i2 (\theta) <= 0.06,
        % 3.1 <= i3 (\shi) <= 3.14
        % 980 <= i4 (\v_own) <= 1200, 
        % 960 <= i5 (\v_in) <= 1200
        % ****NOTE There was a slight mismatch of the ranges of
        % this i5 input for the conference paper, FM2019 "Star-based Reachability of DNNs"
        lb = [1500; -0.06; 3.1; 980; 960];
        ub = [1800; 0.06; 3.14; 1200; 1200];
    elseif property_id == 4        
        % Input Constraints
        % 1500 <= i1(\rho) <= 1800,
        % -0.06 <= i2 (\theta) <= 0.06,
        % (\shi) = 0
        % 1000 <= i4 (\v_own) <= 1200, 
        % 700 <= i5 (\v_in) <= 800
        lb = [1500; -0.06; 0; 1000; 700];
        ub = [1800; 0.06; 0; 1200; 800];
    else
        error('Invalid property ID');
    end
end

% unsafe region
function U = getSpecs(property_id)

    if property_id == 1 
        % output: [x1 = COC; x2 = Weak Left; x3 = Weak Right; x4 = Strong Left; x5 = Strong Right]
        % verify safety: COC <= 1500 or x1 <= 1500 after normalization
        
        % safe region before normalization
        % x1' <= (1500 - 7.5189)/373.9499 = 3.9911
        U = HalfSpace([-1 0 0 0 0], -3.9911); % unsafe region x1' > 3.9911
        
    elseif property_id == 2
        
        % output: [x1 = COC; x2 = Weak Left; x3 = Weak Right; x4 = Strong Left; x5 = Strong Right]
        % safety property: COC is not the maximal score
        % unsafe region: COC is the maximal score: x1 >= x2; x1 >= x3; x1 >= x4, x1
        % >= x5
        unsafe_mat = [-1 1 0 0 0; -1 0 1 0 0; -1 0 0 1 0; -1 0 0 0 1];
        unsafe_vec = [0; 0; 0; 0];
        U = HalfSpace(unsafe_mat, unsafe_vec);
        
    elseif property_id == 3 || property_id == 4
        
        % output: [x1 = COC; x2 = Weak Left; x3 = Weak Right; x4 = Strong Left; x5 = Strong Right]
        % safety property: COC is not the minimal score
        % unsafe region: COC is the minimal score: x1 <= x2; x1 <= x3; x1 <= x4, x1
        % <= x5
        unsafe_mat = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
        unsafe_vec = [0; 0; 0; 0];
        U = HalfSpace(unsafe_mat, unsafe_vec);
        
    else
        error('Invalid property ID');
    end
end