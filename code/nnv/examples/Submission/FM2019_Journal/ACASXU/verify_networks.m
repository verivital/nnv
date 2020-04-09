% verify safety of all 45 ACAS Xu DNNs  

%% verify P3
N1 = 5;
N2 = 9;
[Res, VT] = verify(N1,N2,3,'exact-star',6);
save P3_exact_star.mat Res VT;
[Res, VT] = verify(N1,N2,3,'approx-star',1);
save P3_approx_star.mat Res VT;
[Res, VT] = verify(N1,N2,3, 'approx-zono',1);
save P3_approx_zono.mat Res VT;
[Res, VT] = verify(N1,N2,3, 'abs-dom',1);
save P3_abs_dom.mat Res VT;
% 
%% verify P4
N1 = 5;
N2 = 9;
[Res, VT] = verify(N1,N2,4,'exact-star',6);
save P4_exact_star.mat Res VT;
[Res, VT] = verify(N1,N2,4,'approx-star',1);
save P4_approx_star.mat Res VT;
[Res, VT] = verify(N1,N2,4, 'approx-zono',1);
save P4_approx_zono.mat Res VT;
[Res, VT] = verify(N1,N2,4, 'abs-dom',1);
save P4_abs_dom.mat Res VT;
%% Generate table
generate_tex_P3;
generate_tex_P4;

%% functions
function [Res, VT] = verify(N1, N2, property_id, method, numCores)
s1 = sprintf("\n==================VERIFY PROPERTY P%d==================", property_id);
s2 = sprintf('P%d-%s.txt',property_id, method);
fid = fopen(s2, 'wt');
fprintf(s1);
safes = cell(N1, N2);
vts = zeros(N1, N2);
total_vt = 0;
n_safe = 0;
n_unsafe = 0;
for i=1:N1
    for j=1:N2
        fprintf("\n====================Network N%d_%d=====================", i, j);
        [safe, vt] = verifyNetwork(i, j, property_id, method, numCores);
        if safe == 1
            safes{i, j} = 'UNSAT';
            n_safe = n_safe + 1;
        elseif safe == 0
            safes{i, j} = 'SAT';
            n_unsafe = n_unsafe + 1;
        elseif safe == 2
            safes{i, j} = 'UNK';
        end
        vts(i, j) = vt;
        total_vt = vt + total_vt;
    end
end

% print results
fprintf('\n=====VERIFICATION RESULTS FOR PROPERTY P%d=====', property_id);
fprintf(fid,'\n=====VERIFICATION RESULTS FOR PROPERTY P%d=====', property_id);
fprintf('\n Reachability method: %s', method);
fprintf(fid, '\n Reachability method: %s', method);
fprintf('\n number of cores: %d', numCores);

fprintf('\n Network     Safety     Verification Time');
fprintf(fid,'\n Network     Safety     Verification Time');
for i=1:N1
    for j=1:N2
        fprintf('\n  N%d_%d        %s          %.3f', i, j, safes{i, j}, vts(i, j));
        fprintf(fid, '\n  N%d_%d        %s          %.3f', i, j, safes{i, j}, vts(i, j));
    end
end
fprintf('\nTotal verification time: %.3f', total_vt);
fprintf(fid, '\nTotal verification time: %.3f', total_vt);
fprintf('\nNumber of UNSAT: %d/%d', n_safe, N1*N2);
fprintf(fid,'\nNumber of UNSAT: %d/%d', n_safe, N1*N2);
fprintf('\nNumber of SAT: %d/%d', n_unsafe, N1*N2);
fprintf(fid, '\nNumber of SAT: %d/%d', n_unsafe, N1*N2);
fprintf('\nNumber of UNKNOWN: %d/%d', N1*N2 - n_safe - n_unsafe, N1*N2);
fprintf(fid,'\nNumber of unknown: %d/%d', N1*N2 - n_safe - n_unsafe, N1*N2);
fprintf('\n=======================END=======================');
fprintf(fid,'\n=======================END=======================');

Res = safes;
VT = vts;

end

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