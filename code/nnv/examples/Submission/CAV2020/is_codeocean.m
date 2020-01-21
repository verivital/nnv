function out = is_codeocean()
    if isfolder('/codeocean-true')
        out = 1;
    else
        out = 0;
    end
end

