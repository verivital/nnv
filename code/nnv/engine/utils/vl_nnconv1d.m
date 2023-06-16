function [y] = vl_nnconv1d(input, Weights, Bias, Stride, PaddingSize, DilationFactor)
% the vl_conv1d function provides the 1d convolution calculation of an
% input with given
%   Weights
%   Bias
%   Stride
%   PaddingSize
%   DilationFactor

%   Author: Neelanjana Pal
%   Date: 01/17/2023
    
    %   effext of dilationFactor on filter
    if DilationFactor ~= 1
        effectiveFilterSize = (size(Weights,1)-1)*DilationFactor + 1;
        newFilter = zeros(effectiveFilterSize,size(Weights,2),size(Weights,3));
        j = 1;
        for i=1:size(Weights,1)
            newFilter(j,:,:) = Weights(i);
            j = j + DilationFactor;
        end
    else
        newFilter = Weights;
    end
     n = length(size(input));
    
    if n == 2
        % effect of padding on input
        newInput = zeros(size(input',1)+sum(PaddingSize),size(input',2));
        newInput(PaddingSize(1)+1:PaddingSize(1)+size(input',1),:) = input';
    
        y = zeros((size(newInput,1)-size(newFilter,1))/Stride+1,size(newFilter,3))';
        temp = zeros(size(newInput,1)-size(newFilter,1)+1,size(newFilter,3))';
        
        % 1d convolution
        for i=1: size(newFilter,3)
                for l = 1 : size(newFilter,2)
                    vecA = conv(newInput(:,l)',flip(newFilter(:,l,i)'),'valid');
                    temp(i,:)= temp(i,:)+ vecA;
                end
            if ~isempty(Bias) && ~ isscalar(Bias)
                temp(:,i) = temp(:,i) + Bias(1,i);
            end
        end
    elseif n==3
        % effect of padding on input
        newInput = zeros(size(input,1),sum(PaddingSize)+size(input,2),size(input,3));
        newInput(:,PaddingSize(1)+1:PaddingSize(1)+size(input,2),:) = input;   
        y = zeros(size(newFilter,3),size(newInput,2)-size(newFilter,1)+1,size(newInput,3)); 

        % 1d convolution
        for i=1: size(input,3)
            for k = 1: size(newFilter,3)
                for l = 1: size(newInput,1)
                    vecA = conv(newInput(l,:,i),flip(newFilter(:,l,k)),'valid');
                    y(k,:,i) = y(k,:,i) + vecA;
                end
            end
        end
        if Bias ~=[]
            for i = 1:size(newFilter,3)
                y(i,:,:) = y(i,:,:)+Bias(i);
            end
        end
    end


end