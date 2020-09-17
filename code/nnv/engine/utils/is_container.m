function out = is_container()
    if isfile('/proc/self/cgroup')
        filetext = fileread('/proc/self/cgroup');
        if contains(filetext,'docker')
            out = 1;
        else
            out = 0;
        end
        % 'docker container detected'
    else
        out = 0;
        % 'docker container not detected'
    end
end

