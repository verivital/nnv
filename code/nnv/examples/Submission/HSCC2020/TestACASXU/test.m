    
    str = sprintf('ACASXU_run2a_1_2_batch_2000.mat');
    load(str);
    Layers = [];
    n = length(b);
    for k=1:n - 1
        bk = cell2mat(b(k));
        Wk = cell2mat(W(k));
        Lk = LayerS(Wk, bk, 'poslin');
        Layers = [Layers Lk];
    end
    bn = cell2mat(b(n));
    Wn = cell2mat(W(n));
    Ln = LayerS(Wn, bn, 'purelin');

    Layers = [Layers Ln];
    F = FFNNS(Layers);

    % Input Constraints
    % 55947.69 <= i1(\rho) <= 60760,
    % -3.14 <= i2 (\theta) <= 3.14,
    %-3.14 <= i3 (\shi) <= -3.14
    % 1145 <= i4 (\v_own) <= 1200, 
    % 0 <= i5 (\v_in) <= 60

    lb = [55947.69; -3.14; -3.14; 1145; 0];
    ub = [60760; 3.14; 3.14; 1200; 60];

    % normalize input
    for k=1:5
        lb(k) = (lb(k) - means_for_scaling(k))/range_for_scaling(k);
        ub(k) = (ub(k) - means_for_scaling(k))/range_for_scaling(k);   
    end

    % unsafe: COC is the maximal score
    unsafe_mat = [-1 1 0 0 0; -1 0 1 0 0; -1 0 0 1 0; -1 0 0 0 1];
    unsafe_vec = [0; 0; 0; 0];
    % unsafe region before scaling
%     unsafe_mat = [-1 0 0 0 0];
%     unsafe_vec = [-3.9911];
    
    B = Box(lb, ub);

    U = HalfSpace(unsafe_mat, unsafe_vec); %unsafe region
    
    [safe, VT, counterInput] = F.verify_MSG2(B, U)

    
    
  