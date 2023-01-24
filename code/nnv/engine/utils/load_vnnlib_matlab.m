function [lb_input, ub_input, output] = load_vnnlib_matlab(propertyFile)
    fileID = fopen(propertyFile,'r'); % open vnnlib file
    tline = fgetl(fileID); % Reading file line-by-line
    % initialize variables to track which variables to update (input or output)
    phase = "start"; % there are four phases in each file (declare inputs, declare outputs, define inputs, define outputs)
    % We know the size of the input from the neural network architecture,
    % but we can just avoid memory allocation and dynamically update the
    % input bounds, reshape afterwards
    property = {};
    output = {}; % define output property in a cell matrix (undefined how we want to express this moving forward)
    property_count = 1;
    keyword = [];
    j = 1;
    while ischar(tline)
%         disp(tline)
        if isempty(tline) || startsWith(tline, ';')    % line contains no information
            % Go to next line (no matter which phase we are in)
            tline = fgetl(fileID);
            continue;
        elseif count(tline, '(') ~= count(tline, ')')
                tline = merge_lines(tline, fileID);
                continue;
        elseif phase == "DeclareInput" % start updating input variables
            % Get input dimensions
            if contains(tline, "declare-const") && contains(tline, "X_")
                dim = dim + 1; % only have seen inputs defined as vectors, so this should work
                % the more general approach would require some extra work, but should be easy as well
            elseif contains(tline, "declare-const") && contains(tline, "Y_")
                lb_input = zeros(dim,1);
                ub_input = zeros(dim,1);
                dim = 0; % reset dimension counter
                phase = "DeclareOutput";
                continue;  % redo this line in correct phase
            end
        elseif phase == "DeclareOutput"
            % get output dimensions
            if contains(tline, "declare-const") && contains(tline, "Y_")
                dim = dim + 1; % only have seen inputs defined as vectors, so this should work
                % the more general approach would require some extra work, but should be easy as well
            elseif contains(tline, "assert")
                lb_output = zeros(dim,1);
                ub_output = zeros(dim,1);
                dim = 1; % reset dimension counter
                phase = "DefineInput";
                continue;  % redo this line in correct phase
            end
        elseif phase == "DefineInput" % This only works for properties with one input set
            % assign values to each input dimension
            if contains(tline, "assert") && (contains(tline, 'and') || contains(tline, 'or')) && contains(tline, 'X_')
                error("Currently do not support multiple input options");
            end
            if contains(tline, ">=") || contains(tline, "<=")
                if contains(tline, "X_")
                    s = split(tline, '(');
                    s = s(3:end);
                    for k=1:length(s)
                        t = split(s{k});
                        dim = split(t{2},'_');
                        dim = str2double(dim{2})+1;
                        value = split(t{3},')');
                        value = str2double(value{1});
                        if contains(t{1},">=")
                            lb_input(dim) = value;
                        else
                            ub_input(dim) = value;
                        end
                    end
                    % assign bounds to each input variable
                else
                    phase = "DefineOutput";
                    continue; % redo this line in correct phase
                end
            end
        elseif phase == "DefineOutput"
            % assign output conditions
            if contains(tline, 'assert')
                if ~isempty(property)
                    output{property_count} = property; % add last property to list of property/conditions to verify
                    property_count = property_count + 1;
                end
%                 property = {}; % cell array of conditions to meet
            end
            if contains(tline, '>=') || contains(tline, '<=')
                property = process_condition(tline);
            end
        else
            % initializing (no phase)
            if contains(tline, "declare-const") && contains(tline, "X_")
                phase = "DeclareInput";
                dim = 0;
                continue; % redo this line in correct phase
            end
        end
        % Go to next line
        tline = fgetl(fileID);
%         disp(tline);
    end % end while loop
    output{property_count} = property;
    fclose(fileID); % close vnnlib file
end % end function

%% Helper Functions

function tline = merge_lines(tline, fileID)
    if count(tline, '(') ~= count(tline, ')')
        nextLine = fgetl(fileID);
        tline = [tline nextLine];
        tline = merge_lines(tline, fileID);
    end
end

function assertion = process_condition(tline)
    brackets = 0; % count open brackets
    s = split(tline);
    keywords = containers.Map('KeyType','double','ValueType','any');
    assertion = string();
    assertion_ = string();
    i = 1;
    while i <= length(s)
        if strcmp(s{i}, '(assert')
            brackets = brackets + 1;
            assertion = assertion + "(";
            assertion_ = assertion_ + "(";
        elseif strcmp(s{i}, '(or')
            if contains(assertion, '>=') || contains(assertion, '<=')
                assertion = assertion + " " + keywords(brackets) + " ";
                assertion_ = assertion_ + " " + keywords(brackets) + " ";
            end
            assertion = assertion + "(";
            assertion_ = assertion_ + "(";
            brackets = brackets + 1;
            keywords(brackets) = "||";
        elseif strcmp(s{i}, '(and')
            if contains(assertion, '>=') || contains(assertion, '<=')
                assertion = assertion + " " + keywords(brackets) + " ";
                assertion_ = assertion_ + " " + keywords(brackets) + " ";
            end
            assertion = assertion + "(";
            assertion_ = assertion_ + "(";
            brackets = brackets + 1;
            keywords(brackets) = "&&";
        elseif contains(s{i}, '>=') || contains(s{i}, '<=')
            comp = s{i};
            temp = string();
            temp_ = string();
            if comp(1) == '('
                comp = comp(2:end); % remove parenthesis
                brackets = brackets+1;
                temp = "(";
                temp_ = "(";
            end
            var1 = s{i+1};
            var2 = s{i+2};
            i = i+2;
            var1 = split(var1, "_");
            var2_all = split(var2, ")");
            var2 = var2_all{1};
            var2_all = var2_all(2:end);
            var2 = split(var2, "_");
            if length(var1) < 2
                error('Unkowns output variable definition')
            end
            if strcmp(comp, '<=')
                temp1 = "ub(";
                temp_1 = "lb(";
            else
                temp1 = "lb(";
                temp_1 = "ub(";
            end
            if length(var2) == 1
                temp = temp + temp1 + string(str2double(var1{2})+1) + ")" + " " + string(comp) + " " + string(var2{1}) + ") ";
                if strcmp(comp, "<=")
                    temp_ = temp_ + temp_1 + string(str2double(var1{2})+1) + ")" + " " +  ">=" + " " + string(var2{1}) + ") ";
                else
                    temp_ = temp_ + temp_1 + string(str2double(var1{2})+1) + ")" + " " + "<=" + " " + string(var2{1}) + ") ";
                end
            else
                if strcmp(comp, '>=')
                    temp2 = "ub(";
                    temp_2 = "lb(";
                else
                    temp2 = "lb(";
                    temp_2 = "ub(";
                end
                temp = temp + temp1 + string(str2double(var1{2})+1) + ")" + " " + string(comp) + " " + temp2 + string(str2double(var2{2})+1) + ") ";
                if strcmp(comp, "<=")
                    temp_ = temp_ + temp_1 + string(str2double(var1{2})+1) + ")" + " " + ">="  + " " + temp_2 + string(str2double(var2{2})+1) + ") ";
                else
                    temp_ = temp_ + temp_1 + string(str2double(var1{2})+1) + ")" + " " + "<="  + " " + temp_2 + string(str2double(var2{2})+1) + ") ";
                end
            end
            assertion = assertion + temp;
            assertion_ = assertion_ + temp_;
            par_diff = count(assertion, "(") - count(assertion, ")"); % number of parenthesis remain open
            for k=1:par_diff % add closing brackets to assertion
                assertion = assertion + ")";
                assertion_ = assertion_ + ")";
            end
            brackets = brackets - par_diff; % substract number of closing parenthesis
        end
        i = i+1;
    end
%     disp(assertion);
    assertion = [assertion; assertion_];
end

