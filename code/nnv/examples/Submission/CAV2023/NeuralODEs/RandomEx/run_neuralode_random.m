function run_neuralode_random()

    % Run all random experiments
    
    disp("Extra Small (XS)")
    reach_XS;

    disp("Small (S)")
    reach_S;

    disp("Medium (M)")
    reach_M;

    disp("Large (L)")
    reach_L;

    disp("Extra Large (XL)")
    reach_XL;

    disp("Extra-Extra Large (XXL)")
    reach_XXL;

    disp("Create table results")
    randomEx_table;

end