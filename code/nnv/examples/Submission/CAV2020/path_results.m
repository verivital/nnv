function out = path_results()
    if is_codeocean()
        out = path_results_codeocean();
    else
        out = path_results_normal();
    end
end