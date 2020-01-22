function out = is_codeocean()
    if isfolder('/codeocean-true')
        out = 1;
        % 'codeocean detected'
    else
        out = 0;
        % 'codeocean not detected'
    end
end

