function fcn = layer_fcn(activation_fcn)
    fcn = activation_fcn;
    if fcn == 'relu  '
        fcn = 'ReLU';
    elseif fcn == 'linear'
        fcn = 'Linear';
    end
end

