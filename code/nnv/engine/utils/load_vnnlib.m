function spec = load_vnnlib(propertyFile)

    fileID = fopen(propertyFile,'r'); % open vnnlib file
    tline = fgetl(fileID); % Reading file line-by-line
    phase = "start"; % there are four phases in each file (declare inputs, declare outputs, define inputs, define outputs)
    spec = struct; % Specification, returned as a struct, could also be a list of structs (fields: lb, ub, halfspace)
    property_count = 1;
    
    while ischar(tline)
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
                    property.lb = lb_input;
                    property.ub = ub_input;
                    property.Hg = []; % initialize output of property (to be halfspace)
                    H = []; % half-space matrix
                    g = []; % half-space vector
                    continue; % redo this line in correct phase
                end
            end
        elseif phase == "DefineOutput"
            % assign output conditions
            if contains(tline, 'assert')
                if ~isempty(property)
%                     output{property_count} = property; % add last property to list of property/conditions to verify
                    property_count = property_count + 1;
                end
%                 property = {}; % cell array of conditions to meet
            end
            if contains(tline, '>=') || contains(tline, '<=')
                [H,g] = process_output(tline);
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
    end % end while loop
    spec = property;
    fclose(fileID); % close vnnlib file
end % end function

%% Helper Functions

% Notes:
% Need to add a combination phase for properties that define inputs and outputs together (multiple properties within a property)

function tline = merge_lines(tline, fileID)
    if count(tline, '(') ~= count(tline, ')')
        nextLine = fgetl(fileID);
        tline = [tline nextLine];
        tline = merge_lines(tline, fileID);
    end
end

function [H,g] = process_constraint(tline, dim)
    
end

function [H,g] = process_argument(tline, H, g)
    % if argument contains Y, then add index to H, otherwise add value to g
    
end

function [H,g] = process_condition(tline, dim)
    % parse one assert statement (single line or multiple, but is has been
    % transformed into a single line regardless)
    % tline = char type, text description of assertion to parse
    % dim = output dimensions (type int)
    
    if startsWith(tline,'(<=') || startsWith(tline,'(>=')
       % get the mat and vector of one constraint
       % if a single constraint, H would be a [1,dim] mat, and g a scalar
        [H,g] = process_constraint(tline, dim);
        
    elseif startsWith(tline,'(or')
        % parse all constraints with the "or" statement
        % if multiple *or* statements, H -> cell array of [1,dim], g -> cell array of scalar

    elseif startsWith(tline,'(and')
        % parse all constraints with the "and" statement
        % if multiple *and* statements, then H -> [# ands, dim], and g a [# ands, 1]

    end
end

