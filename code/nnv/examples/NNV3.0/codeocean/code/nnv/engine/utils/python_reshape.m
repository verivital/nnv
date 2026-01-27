function out = python_reshape(input, targetDim)
    
    if length(targetDim) == 3
        temp = reshape(input,targetDim);
        for i = 1 : targetDim(3)
            out(:,:,i) = temp(:,:,i)';
        end
    elseif length(targetDim) == 4 && size(input,2) > 1
        for i = 1 : targetDim(4)
            temp(:,:,:,i) = reshape(input(:,i),targetDim);
            for j = 1 : size(temp(:,:,:,i),3)
                out(:,:,j,i) = temp(:,:,j,i)';
            end
        end
    end

end

