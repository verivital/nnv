function property = load_vnnlib(propertyFile)
    % property = load_vnnlib(propertyFile)
    % input: vnnlib file describing property to verify
    % output:
    %     - property (struct) with following fields
    %         - lb    -> input lower bound vector
    %         - ub    -> input upper bound vector
    %         - prop  -> collection of HalfSpaces describing output specification to verify

%   Assumptions for this function
%     1) Unique input set, only one value for the upper and lower bounds of the inputs
%     2) Output property
%         a) Series of assertions with no "and", "or" conditions                ---->  property.prop: (1x1) HalfSpace
%         b) One assertion with one or more "and", no "or" included             ---->  property.prop: (1x1) HalfSpace
%         c) One assertion with 1 "or" statement followed by N "and" statements ---->  property.prop: (Nx1) HalfSpace
%     3) Constraints must be linear, only >= or <= are supported
%
%    If the vnnlib does not meet these conditions, this function may not work

    fileID = fopen(propertyFile,'r'); % open vnnlib file
    tline = fgetl(fileID); % Reading file line-by-line
    phase = "start"; % there are four phases in each file (declare inputs, declare outputs, define inputs, define outputs)
    
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
                lb_template = zeros(dim,1,'single');
                ub_template = zeros(dim,1,'single');
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
                output_dim = dim;
                dim = 1; % reset dimension counter
                phase = "DefineInput";
                % Initialize variables
                lb_input = lb_template;
                ub_input = ub_template;
                continue;  % redo this line in correct phase
            end
        elseif phase == "DefineInput" % assign lower and upper bounds to input variables (X)
            % There are 3 options (as far as I am aware)
            % 1) One input set -> multiple assert statements (2 per dimension (lower + upper))
            % 2) Multiple input sets, one output set -> or statement defining in same line all constraints for all input dimensions)
            % 3) Multiple input and output sets -> same as above, but including the output constraints in the same line as well

            % What should the output look like?
            % 1) [lb,ub] -> cell array of a 1 vector for lb and ub (e.g.lb = 1x1 cell)
            % 2) lb = cell array of 1xM, M = number of statements in "or"
            % 3) lb = ub = cell array of 1xM, prop = 1xM cell array of props

            % In the case that only some variables can take on multiple
            % input values, we'll throw an error for now (have not seen this case yet)
            
            % Determine which options to go with (option 1 is most common)
            if contains(tline, "assert") && contains(tline, 'or') && contains(tline, 'X_')
                if contains(tline, 'Y_') % option 3
                    [lb_array, ub_array, prop_array] = process_combined_input_output(tline, lb_input, ub_input, output_dim);
                    property.lb = lb_array;
                    property.ub = ub_array;
                    property.prop = prop_array;
                else % option 2
                    [lb_array, ub_array] = process_multiple_inputs(tline, lb_input, ub_input); 
                end
            end

            if contains(tline, ">") || contains(tline, "<") 
                if contains(tline, "X_") % option 1
                    [lb_input, ub_input] = process_input_constraint(tline, lb_input, ub_input);
                   
                else % move onto the next phase (define output constraints)
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
                % There are a few options here
                % 1) multiple single assertions -> One halfspace with M rows (H and g)
                % 2) one or statement with multiple and assertions -> one halfspace for each "and" statement
                %    2a) single assertions -> M halfsapces with only one ow each (H, g)
                %    2b) multiple assertion -> treat each "and" statement as option 1
                
                % check output options

                % Process any output assertion
                if contains(tline, '>=') || contains(tline, '<=') || contains(tline, '>') || contains(tline, '<')
                    ast = process_assertion(tline, output_dim);
                else
                    error("Property not supported yet for assertion: " + string(tline));
                end
                % add assertion (ast) to the property
                if isempty(property.prop)
                    property.prop = {ast};
                else
                    [~, n2] = size(property.prop); % one or multiple halfspace
                    last_ast = property.prop{n2};
                    if length(last_ast) > 1 % means previous ast was an "or"
                        property.prop{n2+1} = ast;
                    else
                        last_ast.Hg.G = [last_ast.Hg.G; ast.Hg.G];
                        last_ast.Hg.g = [last_ast.Hg.g; ast.Hg.g];
                        property.prop{n2} = last_ast;
                    end
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
    fclose(fileID); % close vnnlib file
    % Check if there were multiple inputs
    if exist('lb_array','var')
        property.lb = lb_array;
        property.ub = ub_array;
    end
end % end function

%% Helper Functions

% combine multiple lines into a single one if they belong to a single statement/assertion
function tline = merge_lines(tline, fileID)
    if count(tline, '(') ~= count(tline, ')')
        nextLine = fgetl(fileID);
        tline = [tline nextLine];
        tline = merge_lines(tline, fileID);
    end
end


%% Helper Functions - OUTPUT

% process output assertion
function ast = process_assertion(tline, dim)
    ast = struct;           % initialize assertion
    ast.dim = dim;          % output dimension
    ast.Hg = {};            % collection of HalfSpaces
    tline = strtrim(tline); % remove leading and trailing whitespaces
    % start processing assertion
    tline = tline(8:end); % remove "(assert" from tline
    while tline
        if length(tline) == 1 && tline == ')'
            break
        end
        tline = strtrim(tline);
        % process linear constraint
        if startsWith(tline,'(<=') || startsWith(tline,'(>=') || startsWith(tline,'(<') || startsWith(tline,'(>')
            % get the mat and vector of one constraint
            % if a single constraint, H would be a [1,dim] mat, and g a scalar
            % These get added to ast.Hg
            [ast, len] = process_constraint(tline, ast);
            tline = tline(len+1:end); % update tline, process remainder of assertion
            if isempty(ast.Hg)
                ast.Hg = HalfSpace(ast.H, ast.g);
            else
                [n1, n2] = size(ast.Hg); % size of cell array
                if all([n1 n2 1] == 1)
                    % only one halfspace yet, we can concatenate the assertion here (similar to and statement)
                    ast.Hg.G = [ast.Hg.H; ast.H ];
                    ast.Hg.g = [ast.Hg.g; ast.g];
                else
                    error("Error or property not supported after processing a constraint.");
                end

            end % close
        elseif startsWith(tline,'(or')
            % parse all constraints with the "or" statement
            % if multiple *or* statements, H -> cell array of [1,dim], g -> cell array of scalar
            tline = tline(4:end);
%             pars = 1; % keep track of parenthesis opened during "or" statement
            [ast, tline] = process_or(tline, ast); % there is only one allowed, nothing can come after closing the or statement
            
        elseif startsWith(tline,'(and')
            % parse all constraints with the "and" statement
            % if multiple *and* statements, then H -> [# ands, dim], and g a [# ands, 1]
            error('Property not supported for now. Not allowed -> assertion starting with AND and no OR.')

        elseif startsWith(tline, ')')
            tline = tline(2:end);
        else
            error("Not sure what is happenning, but you should not be here");
        end
    end
end

% process individual output constraint
function [ast, len] = process_constraint(tline, ast)
    % process a single constraint
    len = 1;
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
    H = zeros(1, ast.dim);
    g = 0; 
    if contains(op,"<=") || contains(op, "<") % no need to change signs, same as in halfspace
        H(idx1) = 1;
        if contains(var2,'Y')
            var2 = split(var2, ')');
            var2 = var2{1};
            idx2 = split(var2,'_');
            idx2 = str2double(idx2{2})+1;
            H(idx2) = -1;
        else
            var2 = split(var2, ')');
            g = single(str2double(var2{1}));
        end
    else
        H(idx1) = -1;
        if contains(var2,'Y')
            var2 = split(var2, ')');
            var2 = var2{1};
            idx2 = split(var2,'_');
            idx2 = str2double(idx2{2})+1;
            H(idx2) = 1;
        else
            var2 = split(var2, ')');
            g = -single(str2double(var2{1}));
        end
    end
    % Add constraint (H, g) to assertion variable (ast)
    ast.H = H;
    ast.g = g;
end

% process output "and "statement
function [temp_ast, tline] = process_and(tline, ast)
    temp_ast.dim = ast.dim;
    temp_ast.Hg = {};
    params = 1; % from the removed '(and'
    while ~isempty(tline) && params > 0
        tline = strtrim(tline);
        if startsWith(tline,'(<=') || startsWith(tline,'(>=')
            [temp_ast, len] = process_constraint(tline, temp_ast);
            tline = tline(len+1:end); % update tline, process remainder of assertion
            if isempty(temp_ast.Hg)
                temp_ast.Hg = HalfSpace(temp_ast.H, temp_ast.g);
            else
                temp_ast.Hg.G = [temp_ast.Hg.G; temp_ast.H];
                temp_ast.Hg.g = [temp_ast.Hg.g; temp_ast.g];
            end
            
        elseif startsWith(tline,'(or')
            error("Currently we do not support an OR statement within and AND statement. ")
    %             
        elseif startsWith(tline,'(and')
            error("Currently we do not support an AND statement within and AND statement. ")

        elseif startsWith(tline, ')')
            params = params - 1;
        else
            error("We may be doing something wrong while processing the AND statement or the property is currently not supported.")
        end
    end
end

% process output "or" statement
function [ast, tline] = process_or(tline, ast)
    temp_ast.dim = ast.dim;
    temp_ast.Hg = {};
    pars = 1;
    while ~isempty(tline) > 0 && pars > 0
        tline = strtrim(tline);
        if startsWith(tline,'(<=') || startsWith(tline,'(>=')
            % get the mat and vector of one constraint
            % if a single constraint, H would be a [1,dim] mat, and g a scalar
            % These get added to ast.Hg
            [or_ast, len] = process_constraint(tline, temp_ast);
            tline = tline(len+1:end); % update tline, process remainder of assertion
            or_ast.Hg = HalfSpace(or_ast.H, or_ast.g);
            if isempty(temp_ast.Hg)
                temp_ast = or_ast; % first HalfSpace within OR statement
            else
                temp_ast.Hg = [temp_ast.Hg; or_ast.Hg]; % add HalfSpace to current assertion
            end
        elseif startsWith(tline,'(or')
            error("Currently we do not support an OR statement within and OR statement.")
    
        elseif startsWith(tline,'(and')
            % parse all constraints with the "and" statement
            % if multiple *and* statements, then H -> [# ands, dim], and g a [# ands, 1]
            tline = tline(5:end);
            pars = pars + 1;
            [or_ast, tline] = process_and(tline, temp_ast);
            if isempty(temp_ast.Hg)
                temp_ast = or_ast; % first HalfSpace within OR statement
            else
                temp_ast.Hg = [temp_ast.Hg; or_ast.Hg]; % add HalfSpace to current assertion
            end
        elseif startsWith(tline, ')')
            pars = pars - 1;
            tline = tline(2:end);
        else
            error("Something went wrong while processing the OR statement.")
        end

    end
    ast = temp_ast; % only allow 1 OR statement, so it should be the one and only assertion in the vnnlib file
end

%% Helper Functions - INPUT

% process input constraint
function [lb_input, ub_input] =  process_input_constraint(tline, lb_input, ub_input)
    s = split(tline, '(');
    s = s(end);
    for k=1:length(s)
        t = split(s{k});
        dim = split(t{2},'_');
        dim = str2double(dim{2})+1;
        value = split(t{3},')');
        value = single(str2double(value{1}));
        if contains(t{1},">=")
            lb_input(dim) = value;
        else
            ub_input(dim) = value;
        end
    end
end

% Process input assertion with combined input and output (or statement)
function [lb_array, ub_array, prop_array] = process_combined_input_output(tline, lb_input, ub_input, output_dim)
    all_cons = extractBetween(tline, '(and', '))');
    n = length(all_cons); % how many [lb,ub,prop] combos do we have?
    lb_array = cell(n,1);
    ub_array = cell(n,1);
    prop_array = cell(n,1);
    for i = 1:length(all_cons)
        combo = all_cons{i};
        combo = combo + ")";
        single_cons = extractBetween(combo, "(", ")");
        % initialize a halfspace for each combo
        H = [];
        g = [];
        for k=1:length(single_cons)
            x = single_cons{k};
            if contains(x, 'X')
                [lb_input, ub_input] =  process_input_constraint(x, lb_input, ub_input);
            else
                [H, g] =  process_output_combo_constraint(x, H, g, output_dim);
            end
        end
        Hg = HalfSpace(H,g);
        prop = struct;
        prop.Hg = Hg;
        prop.H = H;
        prop.g = g;
        lb_array{i} = lb_input;
        ub_array{i} = ub_input;
        prop_array{i} = prop;
    end
    
end 

function [H, g] =  process_output_combo_constraint(tline, H, g, output_dim)
    % This is not ideal, but looks like these types of properties only
    % compare Y with values, so can use same code as for input, but need to
    % generalize this for the future
    Hvec = zeros(1, output_dim);
    s = split(tline, '(');
    s = s(end);
    for k=1:length(s)
        t = split(s{k});
        dim = split(t{2},'_');
        dim = str2double(dim{2})+1;
        value = split(t{3},')');
        value = single(str2double(value{1}));
        if contains(t{1},">=") || contains(t{1}, ">")
            Hvec(dim) = -1;
            gval = -value;
        else
            Hvec(dim) = 1;
            gval = value;
        end
        H = [H; Hvec];
        g = [g; gval];
    end
end

% Process input assertion withmultiple input sets (or statement)
function [lb_array, ub_array] = process_multiple_inputs(tline, lb_input, ub_input)
    tline = tline(12:end); % remove '(assert (or'
    pars = 1;
    lb_array = {};
    ub_array = {};
    arr_count = 0;
    while ~isempty(tline) > 0 && pars > 0
        tline = strtrim(tline);
        if startsWith(tline,'(<=') || startsWith(tline,'(>=')
            % get the mat and vector of one constraint
            % if a single constraint, H would be a [1,dim] mat, and g a scalar
            % These get added to ast.Hg
            constraint_char = split(tline, ')');
            constraint_char = constraint_char{1};
            [lb_input, ub_input] = process_input_constraint(constraint_char, lb_input, ub_input);
            tline = tline(length(constraint_char)+2:end); % update tline, process remainder of assertion
        
        elseif startsWith(tline,'(or')
            error("Currently we do not support an OR statement within and OR statement.")
    
        elseif startsWith(tline,'(and')
            % parse all constraints with the "and" statement
            % if multiple *and* statements, then H -> [# ands, dim], and g a [# ands, 1]
            arr_count = arr_count + 1;
            tline = tline(5:end);
            pars = pars + 1;
        
        elseif startsWith(tline, ')')
            pars = pars - 1;
            tline = tline(2:end);
            lb_array{arr_count} = lb_input;
            ub_array{arr_count} = ub_input;
        else
            error("Something went wrong while processing vnnlib property with multiple inputs.")
        end
    end
end
