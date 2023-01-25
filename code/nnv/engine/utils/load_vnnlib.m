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
                    property.prop = {ast};
                else
                    [~, n2] = size(property.prop); % one or multiple halfspace
                    last_ast = property.prop{n2};
                    if length(last_ast) > 1 % means previous ast was an "or"
                        property.prop{n2+1} = ast;
                    else
                        if ast.dim == last_ast.dim
                            last_ast.Hg.G = [last_ast.Hg.G; ast.Hg.G];
                            last_ast.Hg.g = [last_ast.Hg.g; ast.Hg.g];
                            property.prop{n2} = last_ast;
                        else
                            property.prop{n2+1} = ast;
                        end
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



%%%%%%%%%%%%%%%%%%%%%%
%%%%  ASSERTION   %%%%
%%%%%%%%%%%%%%%%%%%%%%
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
        if startsWith(tline,'(<=') || startsWith(tline,'(>=')
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


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   CONSTRAINT  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
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
    if contains(op,"<=") % no need to change signs, same as in halfspace
        H(idx1) = 1;
        if contains(var2,'Y')
            var2 = split(var2, ')');
            var2 = var2{1};
            idx2 = split(var2,'_');
            idx2 = str2double(idx2{2})+1;
            H(idx2) = -1;
        else
            var2 = split(var2, ')');
            g = str2double(var2{1});
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
            g = -str2double(var2{1});
        end
    end
    % Add constraint (H, g) to assertion variable (ast)
    ast.H = H;
    ast.g = g;
end



%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    AND    %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
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



%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    OR    %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
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

