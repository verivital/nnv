star_data = [
       "experiments_results/mnist-s-adv0.1__vrf_result.mat", ...
       "experiments_results/mnist-l-adv0.1-cbd3-vrf__vrf_result.mat", ...
       "experiments_results/mnist-s-adv0.3__vrf_result.mat", ...
       "experiments_results/mnist-l-adv0.3-cbd3-vrf__vrf_result.mat", ...
       "experiments_results/cifar10-s-adv2__vrf_result.mat", ...
       "experiments_results/cifar10-l-adv2-cbd3-vrf__vrf_result.mat", ...
       "experiments_results/cifar10-s-adv8__vrf_result.mat", ...
       "experiments_results/cifar10-l-adv8-cbd3-vrf__vrf_result.mat"
    ];

eevbnn_data = [
       "experiments_results/mnist-s-adv0.1.csv", ...
       "experiments_results/mnist-l-adv0.1-cbd3-vrf.csv", ...
       "experiments_results/mnist-s-adv0.3.csv", ...
       "experiments_results/mnist-l-adv0.3-cbd3-vrf.csv", ...
       "experiments_results/cifar10-s-adv2.csv", ...
       "experiments_results/cifar10-l-adv2-cbd3-vrf.csv", ...
       "experiments_results/cifar10-s-adv8.csv", ...
       "experiments_results/cifar10-l-adv8-cbd3-vrf.csv"
    ];

models_num = 4;

template = "\\multirow{2}{*}{\\textbf{MODEL_NAME}} & DISTURBANCE1 & EEVBNN_TIME_UNSAT1 & EEVBNN_TIME_SAT1 & EEVBNN_NUM_UNSAT1 & EEVBNN_NUM_SAT1 & STAR_TIME_UNSAT1 & STAR_TIME_SAT1 & STAR_TIME_UNK1 & STAR_NUM_UNSAT1 & STAR_NUM_SAT1 & STAR_NUM_UNK1 \\\\ \n & DISTURBANCE2 & EEVBNN_TIME_UNSAT2 & EEVBNN_TIME_SAT2 & EEVBNN_NUM_UNSAT2 & EEVBNN_NUM_SAT2 & STAR_TIME_UNSAT2 & STAR_TIME_SAT2 & STAR_TIME_UNK2 & STAR_NUM_UNSAT2 & STAR_NUM_SAT2 & STAR_NUM_UNK2 \\\\";
      

current_counter = 1;
result = "\\begin{table}[]\n\\centering\n\\resizebox{\\textwidth}{!}{\n\\begin{tabular}{l@{\\hskip 0.2in}l@{\\hskip 0.2in}cc@{\\hskip 0.15in}cc@{\\hskip 0.2in}cс@{\\hskip 0.15in}ccc} \\hline\n\\multirow{3}{*}{\\textbf{Network}} & \\multirow{3}{*}{$\\delta$}& \\multicolumn{4}{c}{\\textbf{EEVBNN}} & \\multicolumn{6}{c}{\\textbf{Approx-Star}} \\\\ \n& & \\multicolumn{2}{l}{Time(s)} & \\multicolumn{2}{l}{\\#Sol} & \\multicolumn{3}{c}{Time(s)} & \\multicolumn{3}{c}{\\#Sol} \\\\  & & UN & S & UN & S & UN & S & UK & UN & S & UK \\\\ \\hline\n";

for i = 1:models_num
    
    
    result = strcat(result, template, " \\hline\n ");
    
    current_path = eevbnn_data(current_counter);
    
    eevbnn_first = csvread(eevbnn_data(current_counter));
    starbnn_first = load(star_data(current_counter)).current_result;
    
    eevbnn_second = csvread(eevbnn_data{current_counter + 1});
    starbnn_second = load(star_data{current_counter + 1}).current_result;
    
    result = strrep(result, "MODEL_NAME", extractBetween(current_path, strfind(current_path, "_result") + 9, strfind(current_path, ".csv") - 1));
    
    result = strrep(result, "DISTURBANCE1", num2str(starbnn_first(1, 2)));
        
    result = strrep(result, "EEVBNN_TIME_UNSAT1", num2str( round(sum(eevbnn_first(eevbnn_first(:,1) == 1, 2)) / length(find(eevbnn_first(:,1) == 1)), 3) ));
    result = strrep(result, "EEVBNN_TIME_SAT1", num2str( round(sum(eevbnn_first(eevbnn_first(:,1) == 0, 2)) / length(find(eevbnn_first(:,1) == 0)), 3) ));
    result = strrep(result, "EEVBNN_NUM_UNSAT1", num2str(length(find(eevbnn_first(:,1) == 1))));
    result = strrep(result, "EEVBNN_NUM_SAT1", num2str(length(find(eevbnn_first(:,1) == 0))));
    
    result = strrep(result, "STAR_TIME_UNSAT1", num2str( round(sum(starbnn_first(starbnn_first(:,3) == 1, 4)) / length(find(starbnn_first(:,3) == 1)), 3) ));
    result = strrep(result, "STAR_TIME_SAT1", num2str( round(sum(starbnn_first(starbnn_first(:,3) == 0, 4)) / length(find(starbnn_first(:,3) == 0)), 3) ));
    result = strrep(result, "STAR_TIME_UNK1", num2str( round(sum(starbnn_first(starbnn_first(:,3) == -1, 4)) / length(find(starbnn_first(:,3) == -1)), 3) ));
    result = strrep(result, "STAR_NUM_UNSAT1", num2str(length(find(starbnn_first(:,3) == 1))));
    result = strrep(result, "STAR_NUM_SAT1", num2str(length(find(starbnn_first(:,3) == 0))));
    result = strrep(result, "STAR_NUM_UNK1", num2str(length(find(starbnn_first(:,3) == -1))));
    
    result = strrep(result, "DISTURBANCE2", num2str(starbnn_second(1, 2)));
        
    result = strrep(result, "EEVBNN_TIME_UNSAT2", num2str( round(sum(eevbnn_second(eevbnn_second(:,1) == 1, 2)) / length(find(eevbnn_second(:,1) == 1)), 3) ));
    result = strrep(result, "EEVBNN_TIME_SAT2", num2str( round(sum(eevbnn_second(eevbnn_second(:,1) == 0, 2)) / length(find(eevbnn_second(:,1) == 0)), 3) ));
    result = strrep(result, "EEVBNN_NUM_UNSAT2", num2str(length(find(eevbnn_second(:,1) == 1))));
    result = strrep(result, "EEVBNN_NUM_SAT2", num2str(length(find(eevbnn_second(:,1) == 0))));
    
    result = strrep(result, "STAR_TIME_UNSAT2", num2str( round(sum(starbnn_second(starbnn_second(:,3) == 1, 4)) / length(find(starbnn_second(:,3) == 1)), 3) ));
    result = strrep(result, "STAR_TIME_SAT2", num2str( round(sum(starbnn_second(starbnn_second(:,3) == 0, 4)) / length(find(starbnn_second(:,3) == 0)), 3) ));
    result = strrep(result, "STAR_TIME_UNK2", num2str( round(sum(starbnn_second(starbnn_second(:,3) == -1, 4)) / length(find(starbnn_second(:,3) == -1)), 3) ));
    result = strrep(result, "STAR_NUM_UNSAT2", num2str(length(find(starbnn_second(:,3) == 1))));
    result = strrep(result, "STAR_NUM_SAT2", num2str(length(find(starbnn_second(:,3) == 0))));
    result = strrep(result, "STAR_NUM_UNK2", num2str(length(find(starbnn_second(:,3) == -1))));
    
    current_counter = current_counter + 2;
end

result = strcat(result, "\\end{tabular}} \n \\caption{ Verification results for MLP1-4. \\textit{Notation}:}\\label{tab:eevbnn_star_001} \n \\end{table}");

fid = fopen('tab_EEVBNN_STAR_comparison.txt','wt');
fprintf(fid, result);
fclose(fid);


starbnn_ver_imgs = load("experiments_results/mnist-l-adv0.1-cbd3-vrf__imgs_result.mat").current_imgs;
starbnn_ver_counterex = load("experiments_results/mnist-l-adv0.1-cbd3-vrf__cexs_result.mat").current_counterexs;



for i=1:5
    csvwrite(strcat('mnist-l-adv0.1_image', num2str(i), '.csv'), starbnn_ver_imgs{i});
    csvwrite(strcat('mnist-l-adv0.1_counterex', num2str(i), '.csv'), starbnn_ver_counterex{i});
end


init_img = starbnn_ver_imgs{1};
counterex = starbnn_ver_counterex{1};

k =0 ;