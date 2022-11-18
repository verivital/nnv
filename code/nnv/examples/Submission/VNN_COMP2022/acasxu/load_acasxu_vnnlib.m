function [lb_input, ub_input, output] = load_acasxu_vnnlib(propertyFile)
    fileID = fopen(propertyFile,'r'); % open vnnlib file
    tline = fgetl(fileID); % Reading file line-by-line
    % initialize variables to track which variables to update (input or output)
    phase = "start"; % there are four phases in each file (declare inputs, declare outputs, define inputs, define outputs)
    % We know the size of the input from the neural network architecture,
    % but we can just avoid memory allocation and dynamically update the
    % input bounds, reshape afterwards
%     lb_input = zeros(5,1);
%     ub_input = zeros(5,1); % define this after declaration
    property = {};
    output = {}; % define output property in a cell matrix (undefined how we want to express this moving forward)
    property_count = 1;
    while ischar(tline)
        if isempty(tline) % line contains no information
            % Go to next line (no matter which phase we are in)
            tline = fgetl(fileID);
            continue;
        elseif phase == "DeclareInput" % start updating input variables
            % Get input dimensions
            if contains(tline, "declare-const") && contains(tline, "X_")
                dim = dm + 1; % only have seen inputs defined as vectors, so this should work
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
                dim = dm + 1; % only have seen inputs defined as vectors, so this should work
                % the more general approach would require some extra work, but should be easy as well
            elseif contains(tline, "assert")
                lb_output = zeros(dim,1);
                ub_output = zeros(dim,1);
                dim = 0; % reset dimension counter
                phase = "DeclareOutput";
                continue;  % redo this line in correct phase
            end
        elseif phase == "DefineInput"
            % assign values to each input dimension
            if contains(tline, ">=") || contains(tline, "<=")
                if contains(tline, "X_")
                    % assign bounds to each input variable
                else
                    phase = "DefineOutput";
                    continue; % redo this line in correct phase
                end
            end
        elseif phase == "DefineOutput"
            % assign output conditions
        else
            % initializing (no phase)
            if contains(tline, "declare-const") && contains(tline, "X_")
                phase = "DeclareInput";
                dim = 0;
                continue; % redo this line in correct phase
            end
        end
%             % nothing else to do in this part
%         elseif contains(tline, 'Output constraints') % stop input, start update output
%             phase = "DefineOutput";
%             i = 1; % reset counter
%         elseif input_start % input structure seems to be well defined and fixed (may only change from X 200 to X_20_10, but never seen that so far)
%             s = split(tline);
%             value = split(s{4},')');
%             value = str2double(value{1});
%             if contains(s{2}, '>=') % lower bound assignment
%                 lb_input(i) = value;
%             else
%                 ub_input(i) = value; % upper bound assignment
%                 i = i+1;
%             end
%         elseif output_start  % output is a little more challenging, as multiple conditions may be defined with different operators (<=, >=, OR, AND,...)
%             %TODO (after this is done, this may be used as general parser)
%             % So far, if multiple conditins in same line after assert parenthesis is opened, 
%             % these are connected with AND operator (set at beginning of line). If within 
%             % the same assert parenthesis, different line will link the conditions using OR statements)
%             % Could also see multiple assert conditions (contiguos lines), could treat them as AND, 
%             % or deal with them as multiple separate properties 
%             % (for generlization, I think the latter may be more appropriate)
%             if contains(tline, 'assert')
%                 if ~isempty(property)
%                     output{property_count} = property; % add last property to list of property/conditions to verify
%                 end
%                 property = {}; % cell array of conditions to meet
%             end
%             if contains(tline, '>=') || contains(tline, '<=')
%                 s = split(tline);
% %                 n = length(s);
%                 if contains(s{1},'assert')
%                     s = s(2:end);
%                 elseif isempty(s{1})
%                     s = s(2:end);
%                 end
%                 s = strjoin(s); % expression with just conditions and keywords like "and"
%                 s = split(s, '('); % get number of conditions to evaluate
%                 for ex=1:length(s) % iterate though each expression
%                     if contains(s{ex}, 'or') % have not seen this yet, but just in case throw an error if encountered
%                         error('We do not support the OR operator here, only AND for now.');
%                     elseif contains(s{ex}, 'and')
%                         j = 1; % condition counter for each property
%                     elseif contains(tline, '>=') || contains(tline, '<=')
%                         condition = split(s{ex});
%                         n = length(condition);
%                         if contains(condition{n}, 'Y') % Comparison on the index (first idx value smaller than 2nd)
%                             property{j,1} = 'index';
%                             var1 = condition{n-1};
%                             var1 = split(var1,'_');
%                             property{j,2} = str2double(var1{end});
%                             var2 = split(condition{n},')');
%                             var2 = var2{1};
%                             var2 = split(var2,'_');
%                             property{j,3} = str2double(var2{end});
%                             i = i+1; % move to next property or value within same property
%                         else
%                             property{j,1} = 'value'; % one value of the output must be smaller than a value
%                             var1 = condition{n-1};
%                             var1 = split(var1,'_');
%                             property{j,2} = str2double(var1{end});
%                             var2 = split(condition{n},')');
%                             var2 = var2{1};
%                             var2 = split(var2,'_');
%                             property{j,3} = str2double(var2{end});
%                         end
%                     end
%                 end
        % Go to next line
        tline = fgetl(fileID);
    end % end while loop
    fclose(fileID); % close vnnlib file
end % end function

