function [lb_input, ub_input, lb_output, ub_output] = load_collins_vnnlib(propertyFile)
    fileID = fopen(propertyFile,'r');
    tline = fgetl(fileID);
    output_start = false;
    input_start = false;
    if endsWith(propertyFile, "w40.vnnlib")
        lb_input = zeros(20,40);
        ub_input = zeros(20,40);
    else
        lb_input = zeros(20,20);
        ub_input = zeros(20,20);
    end
    lb_output = 0;
    ub_output = 0;
    while ischar(tline)
        if isempty(tline)
            % Go to next line
            tline = fgetl(fileID);
            continue;
        elseif contains(tline, 'Input constraints')
            input_start = true;
            tline = fgetl(fileID);
            i = 1;
            continue;
        elseif contains(tline, 'Output constraints')
            input_start = false;
            output_start = true;
            tline = fgetl(fileID); % move to 2 lines from Output Constraints
            i = 1; % reset counter
        elseif input_start
            s = split(tline);
            value = split(s{4},')');
            value = str2double(value{1});
            if contains(s{2}, '>=') % lower bound assignment
                lb_input(i) = value;
            else
                ub_input(i) = value; % upper bound assignment
                i = i+1;
            end
        elseif output_start % This is not great, but all collins properties seem to be the same so it should work
            s = split(tline);
            if contains(s{3}, '>=')
                value = split(s{5},')');
                value = str2double(value{1});
                lb_output = value;
                value = split(s{8},')');
                value = str2double(value{1});
                ub_output = value;
            end
            output_start = false;
        end
        % Go to next line
        tline = fgetl(fileID);
%         disp(tline);
    end
    fclose(fileID);
end

