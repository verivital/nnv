function p = fit_poly_to_mem_consumption_of_LP()
    % data consists of columns of C, rows of C, number of elements
    % of C and memory consumed by one linear program using that C.
    data = [0,      0,      0    ,  0  ;
            3985,   11851,  0.81 ,  0  ;
            169256, 256,    3.78 ,  0  ;
            19752,  3958,   4.98 ,  0  ;
            19709,  3829,   5.84 ,  0  ;
            101389, 4025,   6    ,  0  ;
            147853, 423,    6.38 ,  8  ;
            19907,  4423,   6.67 ,  17 ;
            102759, 8135,   15.5 ,  0  ;
            17797,  53287,  19.6 ,  0  ;
            16875,  75217,  23   ,  0  ;
            105288, 15724,  41   ,  4  ;
            % 105288, 15723,  36.8 ,  0  ;
            % 214358, 135307, 830  ,  0  ;
            ];
    original_mem_consump = data(:, 3);
    data(:, 3) = data(:, 3) + data(:, 4);   % bias for fitting to quadratic
    data(:, 3) = data(:, 3)*2^30;
    x = data(:, [1 2]); % input is number of columns and rows of C
    y = data(:, 3); % output is memory consumed per LP
    % p = polyfitn(x,y,1);
    p = polyfitn(x,y,2);
    % [original_mem_consump, polyvaln(p, x)/2^30]   % debug
end