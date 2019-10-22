N = 1;
M = 1;

fid = fopen('P1_results.txt', 'wt');
fprintf(fid, 'Network       approx-zono         approx-star            abs-dom');
fprintf(fid, '\n              Safety | VT         Safety | VT         Safety | VT');
fprintf(fid, '\n---------------------------------------------------------------------------');
VT1_total = 0;
VT2_total = 0;
VT3_total = 0;

for i=1:N
    for j=1:M
        
        %str = sprintf('ACASXU_run2a_%d_%d_batch_2000.mat', i, j);
        str = sprintf('ACASXU_run2a_1_1_batch_2000.mat');
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

        % unsafe region before scaling
        unsafe_mat = [-1 0 0 0 0];
        unsafe_vec = [-3.9911];
        
        B = Box(lb, ub);

        U = HalfSpace(unsafe_mat, unsafe_vec); %unsafe region
        k = 5; 
        sens_lb = 0.2; % 20%
        [safe1, VT1, ~] = F.verify_MSG(B, 'approx-zono', k, sens_lb, U);
        [safe2, VT2, ~] = F.verify_MSG(B, 'approx-star', k, sens_lb, U);
        [safe3, VT3, ~] = F.verify_MSG(B, 'abs-dom', k, sens_lb, U);
        
        VT1_total = VT1_total + VT1;
        VT2_total = VT2_total + VT2;
        VT3_total = VT3_total + VT3;
        
        fprintf(fid, '\nN_%d_%d         %s    %5.3f       %s     %5.3f       %s    %5.3f', i, j, safe1, VT1, safe2, VT2, safe3, VT3);
     
    end
end
fprintf(fid, '\n---------------------------------------------------------------------------');
fprintf(fid, '\nTOTAL VT:             %5.3f                %5.3f              %5.3f', VT1_total, VT2_total, VT3_total);

fclose(fid);