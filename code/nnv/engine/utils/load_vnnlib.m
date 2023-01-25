function spec = load_vnnlib(propertyFile)

    fileID = fopen(propertyFile,'r'); % open vnnlib file
    tline = fgetl(fileID); % Reading file line-by-line
    phase = "start"; % there are four phases in each file (declare inputs, declare outputs, define inputs, define outputs)
    spec = struct; % Specification, returned as a struct, could also be a list of structs (fields: lb, ub, halfspace)
    prop_count = 1;
    
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
%                 lb_output = zeros(dim,1);
%                 ub_output = zeros(dim,1);
                output_dim = dim;
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
                    property.prop = {}; % initialize output of property (to be a collection of halfspaces)
                    continue; % redo this line in correct phase
                end
            end
        elseif phase == "DefineOutput"
            % define output conditions/property
            if contains(tline, 'assert')
                % check type of assertion, only linear constraints for now
                if contains(tline, '>=') || contains(tline, '<=')
                    ast = process_assertion(tline, output_dim);
                else
                    error("Property not supported yet for assertion: " + string(tline));
                end
                % add assertion (ast) to the property
                if isempty(property.prop)
                    property.prop{prop_count} = ast;
                    prop_count = prop_count + 1;
                else
%                     if 
                end
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
%     spec = property;
    fclose(fileID); % close vnnlib file
end % end function

%% Helper Functions

% Notes:
% Need to add a combination phase for properties that define inputs and outputs together (multiple properties within a property)

% combine multiple lines into a single one if they belong to a single statement/assertion
function tline = merge_lines(tline, fileID)
    if count(tline, '(') ~= count(tline, ')')
        nextLine = fgetl(fileID);
        tline = [tline nextLine];
        tline = merge_lines(tline, fileID);
    end
end

function ast = process_assertion(tline, dim)
    ast = struct;           % initialize assertion
    ast.dim = dim;          % output dimension
    ast.Hg = {};            % collection of HalfSpaces
    tline = strtrim(tline); % remove leading and trailing whitespaces
    % start processing assertion
    tline = tline(7:end); % remove "(assert" from tline
    while tline
        if length(tline) == 1 && tline == ')'
            break
        end
        tline = strtrim(tline);
        % process linear constraint
        if startsWith(tline,'(<=') || startsWith(tline,'(>=')
            % get the mat and vector of one constraint
            % if a single constraint, H would be a [1,dim] mat, and g a scalar
            % These get added to ast.Hg
            [ast, len] = process_constraint(tline, ast);
            tline = tline(len:end); % update tline, process remainder of assertion
            if isempty(ast.Hg)
                ast.Hg = HalfSpace(ast.H, ast.g);
            else
                [n1, n2] = size(ast.Hg); % size of cell array
                if all([n1 n2 1] == 1)
                    % only one halfspace yet, we can concatenate the assertion here (similar to and statement)
                    ast.Hg.G = [ast.Hg.H; ast.H ];
                    ast.Hg.g = [ast.Hg.g; ast.g];
                elseif n2 > n1
                    % we have an array of Halfspace (which means an "or" property predecesess this property, create new HalfSpace
                    ast.Hg{n2+1} = HalfSpace(ast.H, ast.g);
                else
                    error("We should not be here, go back and fix this error!! No other option to add H and g to property");
                end

            end % close
        elseif startsWith(tline,'(or')
            % parse all constraints with the "or" statement
            % if multiple *or* statements, H -> cell array of [1,dim], g -> cell array of scalar
            tline = tline(4:end);
%             pars = 1; % keep track of parenthesis opened during "or" statement
            [ast, tline] = process_or(tline, ast);

        elseif startsWith(tline,'(and')
            % parse all constraints with the "and" statement
            % if multiple *and* statements, then H -> [# ands, dim], and g a [# ands, 1]
            tline = tline(5:end);
            [ast, tline] = process_and(tline, ast);

        elseif startsWith(tline, ')')
            tline = tline(2:end);
        else
            error("Not sure what is happenning, but you should not be here");
        end
    end
end

function [ast, len] = process_constraint(tline, ast)
    % process a single constraint
    len = 0;
    while tline(len) ~= ')'
        len  = len + 1;
    end
    const = tline(1:len);
    const = split(const); % split constraint, should be op first, then Y_a then Y_b or a number
    op = const{1};
    var1 = const{2};
    idx1 = split(var1,'_');
    idx1 = str2double(idx1{2})+1;
    var2 = const{3};    
    % initialize values
    H = zeros(1, dim);
    g = 0; 
    if contains(op,"<=") % no need to change signs, same as in halfspace
        H(idx1) = 1;
        if contains(var2,'Y')
            idx2 = split(var1,'_');
            idx2 = str2double(idx2{2})+1;
            H(idx2) = -1;
        else
            var2 = split(var2, ')');
            g = str2double(var2{1});
        end
    else
        H(idx1) = -1;
        if contains(var2,'Y')
            idx2 = split(var1,'_');
            idx2 = str2double(idx2{2})+1;
            H(idx2) = 1;
        else
            var2 = split(var2, ')');
            g = -str2double(var2{1});
        end
    end
    % Add constraint (H, g) to assertion variable (ast)
    ast.H = H;
    ast.g = g;
end

function [ast, tline] = process_and(tline, ast)
    temp_ast.dim = ast.dim;
    temp_ast.Hg = {};
    while tline
        tline = strtrim(tline);
        if startsWith(tline,'(<=') || startsWith(tline,'(>=')
            % get the mat and vector of one constraint
            % if a single constraint, H would be a [1,dim] mat, and g a scalar
            % These get added to ast.Hg
            [ast, len] = process_constraint(tline, temp_ast);
            tline = tline(len:end); % update tline, process remainder of assertion
        elseif startsWith(tline,'(or')
            % parse all constraints with the "or" statement
            % if multiple *or* statements, H -> cell array of [1,dim], g -> cell array of scalar
            tline = tline(4:end);
    %             pars = 1; % keep track of parenthesis opened during "or" statement
            [ast, tline] = process_or(tline, temp_ast);
    
        elseif startsWith(tline,'(and')
            % parse all constraints with the "and" statement
            % if multiple *and* statements, then H -> [# ands, dim], and g a [# ands, 1]
            tline = tline(5:end);
            [ast, tline] = process_and(tline, temp_ast);

        elseif startsWith(tline, ')')
            break
        else
            error("We may be doing something wrong while processing the or statement.")
        end
    end
    % add temp_ast to ast
    % TODO
end

function [ast, tline] = process_or(tline, ast)
    temp_ast.dim = ast.dim;
    temp_ast.Hg = {};
    while tline
        tline = strtrim(tline);
        if startsWith(tline,'(<=') || startsWith(tline,'(>=')
            % get the mat and vector of one constraint
            % if a single constraint, H would be a [1,dim] mat, and g a scalar
            % These get added to ast.Hg
            [ast, len] = process_constraint(tline, temp_ast);
            tline = tline(len:end); % update tline, process remainder of assertion
        elseif startsWith(tline,'(or')
            % parse all constraints with the "or" statement
            % if multiple *or* statements, H -> cell array of [1,dim], g -> cell array of scalar
            tline = tline(4:end);
    %             pars = 1; % keep track of parenthesis opened during "or" statement
            [ast, tline] = process_or(tline, temp_ast);
    
        elseif startsWith(tline,'(and')
            % parse all constraints with the "and" statement
            % if multiple *and* statements, then H -> [# ands, dim], and g a [# ands, 1]
            tline = tline(5:end);
            [ast, tline] = process_and(tline, temp_ast);
            
        elseif startsWith(tline, ')')
            break
        else
            error("We may be doing something wrong while processing the or statement.")
        end
    end
    % add temp_ast to ast
    % TODO
end

